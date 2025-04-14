/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright 2007 University of Washington
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "ns3/log.h"
#include "ns3/ipv4-address.h"
#include "ns3/ipv6-address.h"
#include "ns3/address-utils.h"
#include "ns3/nstime.h"
#include "ns3/inet-socket-address.h"
#include "ns3/inet6-socket-address.h"
#include "ns3/socket.h"
#include "ns3/udp-socket.h"
#include "ns3/simulator.h"
#include "ns3/socket-factory.h"
#include "ns3/packet.h"
#include "ns3/uinteger.h"

#include "udp-aggregator.h"
#include "udp_switchml_header.h"
#include <utility>
#include <cstdlib>
#include <ctime>
#include "ns3/rng-seed-manager.h"
#include "ns3/random-variable-stream.h"
#include "ns3/double.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <sstream>

#define LYJ_SAATP 0
#define LYJ_SAATP_RATE 0

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("UdpAggregatorApplication");

NS_OBJECT_ENSURE_REGISTERED (UdpAggregator);

TypeId
UdpAggregator::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::UdpAggregator")
    .SetParent<Application> ()
    .SetGroupName("Applications")
    .AddConstructor<UdpAggregator> ()
    .AddAttribute ("Port", "Port on which we listen for incoming packets.",
                   UintegerValue (9),
                   MakeUintegerAccessor (&UdpAggregator::m_port),
                   MakeUintegerChecker<uint16_t> ())
    .AddTraceSource ("Rx", "A packet has been received",
                     MakeTraceSourceAccessor (&UdpAggregator::m_rxTrace),
                     "ns3::Packet::TracedCallback")
    .AddTraceSource ("RxWithAddresses", "A packet has been received",
                     MakeTraceSourceAccessor (&UdpAggregator::m_rxTraceWithAddresses),
                     "ns3::Packet::TwoAddressTracedCallback")
    .AddAttribute ("Msize", "aggregator number",
                   UintegerValue (256),
                   MakeUintegerAccessor (&UdpAggregator::pool_size),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("MAX_Host_Number", "host number",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpAggregator::max_host_num),
                   MakeUintegerChecker<uint16_t> ())
    .AddAttribute ("Level", "switch level",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpAggregator::switch_level),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("Aggregateid", "Aggregateid",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpAggregator::aggregate_id),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("RemoteAddress", 
                   "The destination Address of the outbound packets",
                   AddressValue (),
                   MakeAddressAccessor (&UdpAggregator::m_peerAddress),
                   MakeAddressChecker ())
    .AddAttribute ("RemotePort", 
                   "The destination port of the outbound packets",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpAggregator::m_peerPort),
                   MakeUintegerChecker<uint16_t> ())
    .AddAttribute ("OnSaatp", 
                   "The destination port of the outbound packets",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpAggregator::on_saatp),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("ONOFFASYCC", 
                   "onoffasycc",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpAggregator::onoffasycc),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("ONOFFPS", 
                   "onoffps",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpAggregator::onoffps),
                   MakeUintegerChecker<uint32_t> ())       
    .AddAttribute ("ONOFFTIMEWINDOW", 
                    "onofftimewindow",
                    UintegerValue (1),
                    MakeUintegerAccessor (&UdpAggregator::onofftimewindow),
                    MakeUintegerChecker<uint32_t> ())         
      .AddAttribute("TimeWindow", 
                      "Time Limit",
                      DoubleValue(0.1),
                      MakeDoubleAccessor(&UdpAggregator::m_timeWindow),
                    MakeDoubleChecker<double>())          
  ;
  return tid;
}

UdpAggregator::UdpAggregator ()
{
  RngSeedManager::SetSeed(1);
  minspace = 0;
  collision_marking = 0;
  m_timeDataFile = "samples.txt";
  m_windowSize = 100;  // 默认窗口大小
  m_ewmaAlpha = 0.2;   // EWMA平滑因子
  NS_LOG_FUNCTION (this);
}

UdpAggregator::~UdpAggregator()
{
  NS_LOG_FUNCTION (this);
  m_socket = 0;
  m_socket6 = 0;
  if(bitmap!=nullptr){
      delete[] bitmap;
  }
}

void
UdpAggregator::DoDispose (void)
{
  NS_LOG_FUNCTION (this);
  Application::DoDispose ();
}

void 
UdpAggregator::StartApplication (void)
{
  NS_LOG_FUNCTION (this);
  collision_marking = 0;
  initmem(pool_size,max_host_num);
  
  // 写入功能
  // m_timeLogFile.open(m_timeDataFile, std::ios::app);
  
  // 读取样本并计算最优时间窗口
  LoadCachedSamples();
  
  minspace = 0;
  if (m_socket == 0)
    {
      TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
      m_socket = Socket::CreateSocket (GetNode (), tid);
      InetSocketAddress local = InetSocketAddress (Ipv4Address::GetAny (), m_port);
      // std::cout<<"m_port "<<m_port<<std::endl;
      if (m_socket->Bind (local) == -1)
        {
          NS_FATAL_ERROR ("Failed to bind socket");
        }
      if (addressUtils::IsMulticast (m_local))
        {
          Ptr<UdpSocket> udpSocket = DynamicCast<UdpSocket> (m_socket);
          if (udpSocket)
            {
              // equivalent to setsockopt (MCAST_JOIN_GROUP)
              udpSocket->MulticastJoinGroup (0, m_local);
            }
          else
            {
              NS_FATAL_ERROR ("Error: Failed to join multicast group");
            }
        }
    }

  if (m_socket6 == 0)
    {
      TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
      m_socket6 = Socket::CreateSocket (GetNode (), tid);
      Inet6SocketAddress local6 = Inet6SocketAddress (Ipv6Address::GetAny (), m_port);
      if (m_socket6->Bind (local6) == -1)
        {
          NS_FATAL_ERROR ("Failed to bind socket");
        }
      if (addressUtils::IsMulticast (local6))
        {
          Ptr<UdpSocket> udpSocket = DynamicCast<UdpSocket> (m_socket6);
          if (udpSocket)
            {
              // equivalent to setsockopt (MCAST_JOIN_GROUP)
              udpSocket->MulticastJoinGroup (0, local6);
            }
          else
            {
              NS_FATAL_ERROR ("Error: Failed to join multicast group");
            }
        }
    }

if (m_socket_for_up.empty() && switch_level !=0)
    {
      TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
      for(auto it = multiple_address_for_up.begin();  it != multiple_address_for_up.end(); it++){
        m_socket_for_up[it->first]= Socket::CreateSocket (GetNode (), tid);
        if (Ipv4Address::IsMatchingType(multiple_address_for_up[it->first]) == true)
        {
          
          if (m_socket_for_up[it->first]->Bind () == -1)
            {
              NS_FATAL_ERROR ("Failed to bind socket");
            }
          m_socket_for_up[it->first]->Connect (InetSocketAddress (Ipv4Address::ConvertFrom(multiple_address_for_up[it->first]), multiple_port_for_up[it->first]));
        }
        else if (Ipv6Address::IsMatchingType(multiple_address_for_up[it->first]) == true)
          {
            if (m_socket_for_up[it->first]->Bind6 () == -1)
              {
                NS_FATAL_ERROR ("Failed to bind socket");
              }
            m_socket_for_up[it->first]->Connect (Inet6SocketAddress (Ipv6Address::ConvertFrom(multiple_address_for_up[it->first]), multiple_port_for_up[it->first]));
          }
        else if (InetSocketAddress::IsMatchingType (multiple_address_for_up[it->first]) == true)
          {
            if (m_socket_for_up[it->first]->Bind () == -1)
              {
                NS_FATAL_ERROR ("Failed to bind socket");
              }
            m_socket_for_up[it->first]->Connect (multiple_address_for_up[it->first]);
          }
        else if (Inet6SocketAddress::IsMatchingType (multiple_address_for_up[it->first]) == true)
          {
            if (m_socket_for_up[it->first]->Bind6 () == -1)
              {
                NS_FATAL_ERROR ("Failed to bind socket");
              }
            m_socket_for_up[it->first]->Connect (multiple_address_for_up[it->first]);
          }
        else
          {
            NS_ASSERT_MSG (false, "Incompatible address type: " << multiple_address_for_up[it->first]);
          }
        }
    }

  m_socket->SetRecvCallback (MakeCallback (&UdpAggregator::HandleRead, this));
  m_socket6->SetRecvCallback (MakeCallback (&UdpAggregator::HandleRead, this));
  if (switch_level !=0){
    for(auto it = multiple_address_for_up.begin();  it != multiple_address_for_up.end(); it++){
      m_socket_for_up[it->first]->SetRecvCallback (MakeCallback (&UdpAggregator::HandleRead, this));
      NS_LOG_LYJ("init socket");
    }
  }


  //  if (m_socket_for_up.empty() && switch_level !=0)
  //   {
  //     TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  //     for(uint32_t i = 0;  i<multiple_address_for_up.size(); i++){
  //       m_socket_for_up.push_back(Socket::CreateSocket (GetNode (), tid));
  //       if (Ipv4Address::IsMatchingType(multiple_address_for_up[i]) == true)
  //       {
  //         if (m_socket_for_up[i]->Bind () == -1)
  //           {
  //             NS_FATAL_ERROR ("Failed to bind socket");
  //           }
  //         m_socket_for_up[i]->Connect (InetSocketAddress (Ipv4Address::ConvertFrom(multiple_address_for_up[i]), multiple_port_for_up[i]));
  //       }
  //       else if (Ipv6Address::IsMatchingType(multiple_address_for_up[i]) == true)
  //         {
  //           if (m_socket_for_up[i]->Bind6 () == -1)
  //             {
  //               NS_FATAL_ERROR ("Failed to bind socket");
  //             }
  //           m_socket_for_up[i]->Connect (Inet6SocketAddress (Ipv6Address::ConvertFrom(multiple_address_for_up[i]), multiple_port_for_up[i]));
  //         }
  //       else if (InetSocketAddress::IsMatchingType (multiple_address_for_up[i]) == true)
  //         {
  //           if (m_socket_for_up[i]->Bind () == -1)
  //             {
  //               NS_FATAL_ERROR ("Failed to bind socket");
  //             }
  //           m_socket_for_up[i]->Connect (multiple_address_for_up[i]);
  //         }
  //       else if (Inet6SocketAddress::IsMatchingType (multiple_address_for_up[i]) == true)
  //         {
  //           if (m_socket_for_up[i]->Bind6 () == -1)
  //             {
  //               NS_FATAL_ERROR ("Failed to bind socket");
  //             }
  //           m_socket_for_up[i]->Connect (multiple_address_for_up[i]);
  //         }
  //       else
  //         {
  //           NS_ASSERT_MSG (false, "Incompatible address type: " << multiple_address_for_up[i]);
  //         }
  //       }
  //   }
  // m_socket->SetRecvCallback (MakeCallback (&UdpAggregator::HandleRead, this));
  // m_socket6->SetRecvCallback (MakeCallback (&UdpAggregator::HandleRead, this));
  // if (switch_level !=0){
  //   for(uint32_t i = 0;  i<multiple_address_for_up.size(); i++){
  //     m_socket_for_up[i]->SetRecvCallback (MakeCallback (&UdpAggregator::HandleRead, this));
  //     NS_LOG_LYJ("init socket");
  //   }
  // }
  // outputthrought();
  Simulator::Schedule(Seconds(1.0), &UdpAggregator::outputthrought, this);
   
  endkey = 0;
}

void 
UdpAggregator::StopApplication ()
{
  NS_LOG_FUNCTION (this);
  if (m_timeLogFile.is_open()) {
    m_timeLogFile.close();
  }
  if (m_socket != 0) 
    {
      m_socket->Close ();
      m_socket->SetRecvCallback (MakeNullCallback<void, Ptr<Socket> > ());
    }
  if (m_socket6 != 0) 
    {
      m_socket6->Close ();
      m_socket6->SetRecvCallback (MakeNullCallback<void, Ptr<Socket> > ());
    }
}

void 
UdpAggregator::UpdateLambdaEstimate(uint16_t appId, Time arrivalTime)
{
  auto& window = m_arrivalWindows[appId];
  window.push_back(arrivalTime);

  if (window.size() > m_windowSize) {
    window.pop_front();
  }

  if (window.size() < 2) {
    return;
  }

  double totalIntervalSeconds = 0.0;
  for (size_t i = 1; i < window.size(); ++i) {
    double intervalSeconds = (window[i] - window[i-1]).GetSeconds();
    // 添加保护，避免间隔时间为0
    if (intervalSeconds < 1e-9) {  // 小于1纳秒的当作1纳秒处理
      intervalSeconds = 1e-9;
    }
    totalIntervalSeconds += intervalSeconds;
  }
  
  double avgIntervalSeconds = totalIntervalSeconds / (window.size() - 1);
  // 添加保护，避免除零
  if (avgIntervalSeconds < 1e-9) {
    avgIntervalSeconds = 1e-9;
  }

  double currentLambda = 1.0 / avgIntervalSeconds;
  
  // 添加上限保护
  if (currentLambda > 1e9) {  // 限制最大到达率
    currentLambda = 1e9;
  }

  if (m_lambda.find(appId) == m_lambda.end()) {
    m_lambda[appId] = currentLambda;
  } else {
    m_lambda[appId] = m_ewmaAlpha * currentLambda + (1 - m_ewmaAlpha) * m_lambda[appId];
  }

  NS_LOG_INFO("App " << appId << " lambda estimate: " << m_lambda[appId] << " packets/sec");


  // std::cout << "App " << appId 
  //           << " lambda=" << m_lambda[appId] 
  //           << " packets/sec, n=" << n 
  //           << " nodes" << std::endl;
}
// 自适应辛普森积分
template<typename Func>
double integrate(Func f, double a, double b, double eps = 1e-6, int max_depth = 20) {
    using AdaptiveFunc = std::function<double(Func, double, double, double, double, double, double, int)>;
    
    AdaptiveFunc adaptive;
    adaptive = [&](Func func, double a, double b, double eps, double S, double fa, double fb, int depth) {
        double c = (a + b) / 2, h = b - a;
        double fc = func(c);
        double d = (a + c) / 2, e = (c + b) / 2;
        double fd = func(d), fe = func(e);
        double Sleft = (c - a) / 6 * (fa + 4*fd + fc);
        double Sright = (b - c) / 6 * (fc + 4*fe + fb);
        double S2 = Sleft + Sright;

        if (depth <= 0 || std::abs(S2 - S) <= 15*eps)
            return S2 + (S2 - S)/15;
        return adaptive(func, a, c, eps/2, Sleft, fa, fc, depth-1) +
               adaptive(func, c, b, eps/2, Sright, fc, fb, depth-1);
    };

    double fa = f(a), fb = f(b);
    return adaptive(f, a, b, eps, (b - a)/6*(fa + 4*f((a+b)/2) + fb), fa, fb, max_depth);
}
// Erlang B公式计算（递推实现）
double erlangB(int Y, double rho) {
  double B = 1.0;
  for (int k = 1; k <= Y; ++k) {
      B = (rho * B) / (k + rho * B);
  }
  return B;
}
double computeExpectedD(int n, double T, double alpha, double mu) {
  auto F = [&](double x) { return 1 - pow(1 + x / mu, -alpha); };
  auto f = [&](double x) { return (alpha / mu) * pow(1 + x / mu, -(alpha + 1)); };

  // 积分项
  double integral = integrate([&](double x) {
      return x * (n - 1) * pow(F(x), n - 2) * f(x);
  }, 0, T);

  // 第二项
  double term2 = T * (1 - pow(F(T), n - 1));
  return integral + term2;
}
double computeS(double T, double lambda, int Y, int n, double alpha, double mu) {
  double ED = computeExpectedD(n, T, alpha, mu);
  double rho = lambda * ED;
  double B = erlangB(Y, rho);
  double F_T = 1 - pow(1 + T / mu, -alpha);
  return lambda * (1 - B) * (1 + (n - 1) * F_T);
}
double goldenSectionSearch(double lambda, int Y, int n, double alpha, double mu,
  double a, double b, double tol = 1e-5) {
    const double ratio = (sqrt(5) - 1) / 2;
    double c = b - ratio * (b - a);
    double d = a + ratio * (b - a);

    while (abs(c - d) > tol) {
      double sc = computeS(c, lambda, Y, n, alpha, mu);
      double sd = computeS(d, lambda, Y, n, alpha, mu);
      if (sc > sd) {
        b = d;
      } else {
        a = c;
      }
      c = b - ratio * (b - a);
      d = a + ratio * (b - a);
    }
    return (a + b) / 2;
}
void 
UdpAggregator::EstimateAndUpdateT(uint16_t appid, uint32_t key) {
  // 转换时间样本为double格式（秒）
  std::vector<double> samples;
  for (auto& t : m_cachedSamples[std::make_pair(appid, key)]) {
      samples.push_back(t);
  }

  // Lomax参数估计
  auto params = estimateLomaxParameters(samples);
  alpha[appid] = params[0];
  mu[appid] = params[1];

  // 黄金分割搜索找最优T
  double Tmax = samples[0];
  for (auto x : samples) {
    Tmax = std::max(Tmax, x);
  }
  double optimalT = goldenSectionSearch(m_lambda[appid], Y, n,
                                      params[0], params[1], 0, Tmax);
  
  // 更新当前T值（时间窗口）
  m_currentT[appid] = optimalT;
  std::cout << "Alpha: " << params[0] << " Mu: " << params[1] << std::endl;
  std::cout << "Optimal Time Window: " << optimalT << " seconds" << std::endl;
  NS_LOG_INFO("Updated T for app " << appid << " to " << optimalT << "s");
}

std::vector<double> 
UdpAggregator::estimateLomaxParameters(const std::vector<double>& samples) {
  int m = samples.size();
  double mu = 1.0;  // 初始猜测
  double alpha = 1.0;
  const double tolerance = 1e-6;
  const int max_iter = 100;

  for (int iter = 0; iter < max_iter; ++iter) {
      // 更新alpha
      double sum_log = 0.0;
      for (double x : samples) sum_log += log(1 + x / mu);
      double alpha_new = m / sum_log;

      // 牛顿法更新mu
      double mu_new = mu;
      const double target = m / (alpha_new + 1);
      for (int newton_iter = 0; newton_iter < max_iter; ++newton_iter) {
          double sum = 0.0, deriv = 0.0;
          for (double x : samples) {
              double denom = mu_new * (mu_new + x);
              sum += x / denom;
              deriv += -x * (2 * mu_new + x) / (mu_new * mu_new * pow(mu_new + x, 2));
          }
          double error = sum - target;
          if (abs(error) < tolerance) break;
          mu_new -= error / deriv;
      }

      // 打印每次迭代的参数更新情况
      // std::cout << "Iteration " << iter << ": alpha = " << alpha_new << ", mu = " << mu_new << std::endl;

      // 检查收敛
      if (abs(mu_new - mu) < tolerance && abs(alpha_new - alpha) < tolerance) {
          // std::cout << "Converged after " << iter << " iterations." << std::endl;
          return {alpha_new, mu_new};
      }
      mu = mu_new;
      alpha = alpha_new;
  }
  // std::cout << "Did not converge within " << max_iter << " iterations." << std::endl;
  return {alpha, mu};
}
void 
UdpAggregator::HandleRead (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);

  Ptr<Packet> packet;
  Address from;
  Address localAddress;
  while ((packet = socket->RecvFrom (from)))
    {
      
      socket->GetSockName (localAddress);
      m_rxTrace (packet);
      m_rxTraceWithAddresses (packet, from, localAddress);
      SwitchHeader tmp_header;
      packet->RemoveHeader (tmp_header);
      uint32_t recv_key = tmp_header.GetKey();
      uint8_t recv_id = tmp_header.GetHostid();
      uint8_t recv_num = tmp_header.GetHostnum();
      uint16_t isAck = tmp_header.GetACK();
      uint16_t isCollision = tmp_header.GetCollision();
     // uint16_t isEnd = tmp_header.GetEnd();
      uint16_t app_id = tmp_header.GetAppID();
      uint8_t is_bpAggr = tmp_header.GetbpAggr();
      uint32_t syn_urgent = tmp_header.GetDelta();
      uint8_t set_ecn = tmp_header.GetMACK();

      // 只处理非背景流量的包
      if (!is_bpAggr && !isAck) {  // 如果不是背景聚合流量且不是ACK
        m_observedHosts[app_id].insert(recv_id);
        Time now = Simulator::Now();
        n = recv_num;  // 使用包中的实际节点数
        
        // 可以添加验证
        if (m_observedHosts[app_id].size() > recv_num) {
          std::cout << "Warning: App " << app_id 
                    << " observed more hosts (" << m_observedHosts[app_id].size() 
                    << ") than reported in packet (" << (int)recv_num << ")" << std::endl;
        }
        
        UpdateLambdaEstimate(app_id, now);
      }
      // SwmlRouteTag temptag;
      // packet->RemovePacketTag(temptag);


      NS_LOG_INFO("rec: ["<<app_id<<", "<<recv_id<<", "<<recv_key<<"], ack :"<<isAck);
      if (isCollision==1 && isAck!=1){
        std::cout<<"debug udp-aggregator isCollision==1 && isAck!=1"<<std::endl;
        m_count_sent[recv_id]++;
        tmp_header.SetCollision(1);
        tmp_header.SetTotal(1);
        tmp_header.SetMerged(0);
        packet->AddHeader(tmp_header);
        m_socket_for_up[app_id]->Send (packet); //to up worker
        
      }
      else if(isAck==1){
        NS_LOG_INFO("aggregator recieve an ack ["<<app_id<<", "<<recv_key<<"]");
        bool mark_col = isCollision;
        if(isexist(app_id,recv_key)){ //check marking
          int tmpindex = app_and_key_to_bitmap_index[app_id][recv_key];
          mark_col = aggr_collision[tmpindex];
        }
        cleanaggregator(app_id,recv_key);
        broadcast(socket,packet,app_id,recv_key, is_bpAggr, mark_col, syn_urgent, set_ecn);
        // std::cout<<"aggregator recieve isack ["<<app_id<<", "<<recv_key<<"]"<<std::endl;
      }else if(switch_level==0){
        std::cout<<"switch_level==0"<<std::endl;
        if(group[app_id].size() != recv_num){
            group[app_id].insert(from);
        }
        NS_LOG_INFO("server get key:"<<recv_key);

        NS_LOG_INFO("server get key:"<<recv_key);

        NS_LOG_INFO(
          "At time "<< Simulator::Now ().As (Time::S) <<
        "server received key "<< recv_key << " from (" <<
        InetSocketAddress::ConvertFrom (from).GetIpv4 () << " , " <<
        InetSocketAddress::ConvertFrom (from).GetPort ()<<")");
        if (InetSocketAddress::IsMatchingType (from))
            {
            NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " server received " << packet->GetSize () << " bytes from " <<
                        InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
                        InetSocketAddress::ConvertFrom (from).GetPort ());
            }
        else if (Inet6SocketAddress::IsMatchingType (from))
            {
            NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " server received " << packet->GetSize () << " bytes from " <<
                        Inet6SocketAddress::ConvertFrom (from).GetIpv6 () << " port " <<
                        Inet6SocketAddress::ConvertFrom (from).GetPort ());
            }

        //packet->RemoveAllPacketTags ();
        packet->RemoveAllByteTags ();
        uint32_t totaldelta = 0;
        if(aggregate_pkt(app_id, recv_key, recv_id, recv_num, packet, socket, from, is_bpAggr, syn_urgent)){
            // broadcast(socket,packet,app_id,recv_key, totaldelta);
            upload(packet,app_id,recv_key, is_bpAggr, syn_urgent);
        }
      }else{
        if(group[app_id].size() != recv_num){
            group[app_id].insert(from);
        }
        //NS_LOG_LYJ("server get key from "<<recv_id<<" in level1:"<<recv_key);
        if (InetSocketAddress::IsMatchingType (from))
            {
            NS_LOG_INFO ("InetSocketAddress::IsMatchingType (from) At time " << Simulator::Now ().As (Time::S) << " server received " << packet->GetSize () << " bytes from " <<
                        InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
                        InetSocketAddress::ConvertFrom (from).GetPort ());
            }
        else if (Inet6SocketAddress::IsMatchingType (from))
            {
            NS_LOG_INFO ("Inet6SocketAddress::IsMatchingType (from) At time " << Simulator::Now ().As (Time::S) << " server received " << packet->GetSize () << " bytes from " <<
                        Inet6SocketAddress::ConvertFrom (from).GetIpv6 () << " port " <<
                        Inet6SocketAddress::ConvertFrom (from).GetPort ());
            }

        //packet->RemoveAllPacketTags ();
        packet->RemoveAllByteTags ();

        m_count_sent[recv_id]++;

        // UDPecn tag;

        uint32_t totaldelta = 0;

        if(aggregate_pkt(app_id, recv_key, recv_id, recv_num, packet, socket, from, is_bpAggr, syn_urgent)){

          //std::cout<<"totaldelta "<<totaldelta<<std::endl;
          if(!onoffps){
            cleanaggregator(app_id,recv_key);
            broadcast(socket,packet,app_id,recv_key, is_bpAggr, isCollision, syn_urgent, set_ecn);
          }else{
            upload(packet,app_id,recv_key, is_bpAggr, syn_urgent);
            if (onofftimewindow) cleanaggregator(app_id, recv_key);
          }
        }
      }

      if(!isAck) {  // 只在收到非ACK包时记录
        // LogArrivalTime(app_id, recv_key);  // Commented out as we're using cached samples
      }
    }
}

void
UdpAggregator::initmem (uint32_t size, uint16_t host)
{
    int app_num = 256;
    pool_size = size;
    max_host_num = host;
    
    bitmap = new bool*[pool_size];
    interval = new int64_t*[pool_size];
    count_pkt = new uint32_t[pool_size];
    aggr_collision = new bool[pool_size];
    idx_appid = new uint16_t[pool_size];
    start_time_index = new int64_t[pool_size];
    delta_host = new uint32_t*[app_num];
    delta_sum = new uint32_t[app_num];

    for(int i = 0; i < app_num; i++){
      delta_sum[i] = 0;
      delta_host[i] = new uint32_t[max_host_num];
      for(uint32_t j = 0; j < max_host_num; j++){
        delta_host[i][j] = 0;
      }
    }


    for(uint32_t i = 0; i < pool_size; i++){
        bitmap[i] = new bool[max_host_num];
        interval[i] = new int64_t[max_host_num];
        for(uint32_t j = 0; j < max_host_num; j++){
            bitmap[i][j] = false;
            interval[i][j] = 0;
        }
        unused.push_back(i);
        count_pkt[i] = 0;
        aggr_collision[i] = 0;
        idx_appid[i] = 65534;
        start_time_index[i] = -1;
        
    }

    for(uint32_t i = 0; i < max_host_num; i++){
      m_count_sent.push_back(0);
      m_count_key.push_back(0);
    }
}
bool
UdpAggregator::aggregate_pkt(uint16_t appid, uint32_t key, uint8_t hostid, uint8_t hostnum, Ptr<Packet> packet, Ptr<Socket> socket, Address from, uint8_t is_bpAggr, uint32_t syn_urgent)
{
    uint32_t index = 0;
    if (appid == 1){
        NS_LOG_INFO("Aggregating packet for app1: host=" << (int)hostid << " key=" << key);
    }
    
    if (onofftimewindow && hasForwarded(appid, key)) {
        SwitchHeader tmp_header;
        tmp_header.SetKey (key);
        tmp_header.SetHostid(hostid);
        tmp_header.SetHostnum(hostnum);
        tmp_header.SetACK(0);
        tmp_header.SetCollision(0);
        tmp_header.SetAppID(appid);
        tmp_header.SetDelta(syn_urgent);
        tmp_header.SetTotal(1);
        tmp_header.SetMerged(0);
        
        packet->AddHeader(tmp_header);
        m_forwarded+=1;
        m_timeout_forwarded+=1;
        m_socket_for_up[appid]->Send (packet);
        return false;
    }

    if (is_bpAggr){

      SwitchHeader tmp_header;
      tmp_header.SetKey (key); //set the key in header
      tmp_header.SetHostid(hostid);
      tmp_header.SetHostnum(hostnum);
      tmp_header.SetACK(0);
      tmp_header.SetCollision(0);
      tmp_header.SetAppID(appid);
      tmp_header.SetbpAggr(1);
      tmp_header.SetDelta(syn_urgent);
      tmp_header.SetTotal(1);
      tmp_header.SetMerged(0);
      packet->AddHeader(tmp_header);

      m_socket_for_up[appid]->Send (packet);
      return false;
    }

    if(isexist(appid,key)){ //already exist
        index = app_and_key_to_bitmap_index[appid][key];

    }else{
      if (iskicked(appid,key)){
        app_and_key_to_kicked[appid][key]++;
        if (app_and_key_to_kicked[appid][key] == hostnum){
          app_and_key_to_kicked[appid].erase(key); //clean kick table to reduce computation time
        }
        // std::cout<<"already be kicked! and ["<<appid<<", "<<(int)hostid<<", "<<key<<"] forward to ps"<<std::endl;
        SwitchHeader tmp_header;
        tmp_header.SetKey (key); //set the key in header
        tmp_header.SetHostid(hostid);
        tmp_header.SetHostnum(hostnum);
        tmp_header.SetACK(0);
        tmp_header.SetCollision(1);
        tmp_header.SetAppID(appid);
        tmp_header.SetDelta(syn_urgent);
        tmp_header.SetTotal(1);
        tmp_header.SetMerged(0);
        packet->AddHeader(tmp_header);

        m_socket_for_up[appid]->Send (packet);
        return false;
      }
      if(unused.empty()){ // no space
        if (onofftimewindow == 0) {
          collision_marking ++;
          collision_marking = collision_marking > pool_size ? pool_size : collision_marking;
  
          app_and_key_to_kicked[appid][key] = 1; //add this key to kick table
        
          
          // printf("debug\n");
          Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable>();
          uint32_t col_index = uv->GetInteger(0,pool_size)%pool_size; //unused
          // printf("col_index %d\n", col_index);
          // int col_index = 1;
          int looptime = pool_size < 100 ? pool_size : 100;
          for (int k = 0; k < looptime && collision_marking > 0; k++){
            col_index = uv->GetInteger(0,pool_size)%pool_size;
            if (idx_appid[col_index] != appid && aggr_collision[col_index] != 1){
              aggr_collision[col_index] = 1;
              collision_marking --;
            }
          }
        }
        
        
        //NS_LOG_LYJ("no aggregators! and ["<<appid<<", "<<host<<", "<<key<<"] forward to ps");
        // std::cout<<"no aggregators! and ["<<appid<<", "<<hostid<<", "<<key<<"] forward to ps"<<std::endl;
        SwitchHeader tmp_header;
        tmp_header.SetKey (key); //set the key in header
        tmp_header.SetHostid(hostid);
        tmp_header.SetHostnum(hostnum);
        tmp_header.SetACK(0);
        tmp_header.SetCollision(1);
        tmp_header.SetAppID(appid);
        tmp_header.SetDelta(syn_urgent);
        tmp_header.SetTotal(1);
        tmp_header.SetMerged(0);
        packet->AddHeader(tmp_header);


        m_forwarded += 1;
        m_socket_for_up[appid]->Send (packet); //to up worker
        //NS_LOG_LYJ("faile! time "<<Simulator::Now ().As (Time::S)<<"rec: ["<<appid<<", "<<host<<", "<<key<<"] update bitmap in index:("<<index<<","<<host<<")");
        
        
        return false;
      }else{ //have space
            //NS_LOG_LYJ("Rest Aggregator: "<<unused.size());
            index = unused.back();
            unused.pop_back();
            updateindexmap(appid,key,index);
            index = app_and_key_to_bitmap_index[appid][key];
            idx_appid[index] = appid;
            AggregatorUnit& unit = m_units[appid][key];
             unit.mergedPacket = packet->Copy();
             unit.startTime = Simulator::Now();
            //  std::cout << onofftimewindow << "\n";
            //  if (onofftimewindow) 
            //      unit.timer = Simulator::Schedule(Seconds(m_timeWindow),
            //                &UdpAggregator:: ForceFlush, this, appid, key);
             if (onofftimewindow) {
              
              double currentT = m_currentT.count(appid) ? m_currentT[appid] : m_timeWindow;
              // std::cout << "m_timewindow" << currentT << std::endl;
              // unit.timer = Simulator::Schedule(Seconds(currentT),
              //               &UdpAggregator::ForceFlush, this, appid, key);

              unit.timer = Simulator::Schedule(Seconds(currentT),
                            &UdpAggregator::ForceFlush, this, appid, key);
             }
            // key_to_bitmap_index.insert(std::pair<uint32_t, uint32_t>(key,index));
            //get time
            if (unused.size() > minspace) minspace = unused.size();
            UDPecn tag;
    
          
            if(packet->RemovePacketTag (tag)){
             NS_LOG_INFO("get direct pkt with ecn: "<<(uint32_t)tag.GetECN());
            }

      }
    }
    NS_LOG_INFO("get index["<<index<<"], but max is "<<pool_size);
    NS_LOG_INFO("time "<<Simulator::Now ().As (Time::S)<<" rec: ["<<appid<<", "<<hostid<<", "<<key<<"] update bitmap in index:("<<index<<","<<hostid<<")");
    //std::cout<<"time "<<Simulator::Now ().As (Time::S)<<" rec: ["<<appid<<", "<<hostid<<", "<<key<<"] update bitmap in index:("<<index<<","<<hostid<<")"<<std::endl;
    if(bitmap[index][hostid]==false){
        bitmap[index][hostid] = true;
        count_pkt[index]++;

        if (appid == 1) {
            NS_LOG_INFO("App1 progress: " << count_pkt[index] << "/" << (int)hostnum);
            
            // 检查是否长时间没有进展
            static Time lastProgressTime = Simulator::Now();
            static uint32_t lastProgressCount = 0;
            
            if (count_pkt[index] > lastProgressCount) {
                lastProgressCount = count_pkt[index];
                lastProgressTime = Simulator::Now();
            } else if ((Simulator::Now() - lastProgressTime).GetSeconds() > 1.0) {
                NS_LOG_WARN("App1 aggregation stalled at " << count_pkt[index] << "/" << (int)hostnum);
                // 强制完成聚合
                if (onofftimewindow) {
                    if (m_units[appid].count(key)) {
                        Simulator::Cancel(m_units[appid][key].timer);
                        m_units[appid].erase(key);
                    }
                }
                return true;
            }
        }

        if(count_pkt[index]== hostnum){
            if (onofftimewindow) {
                if (m_units[appid].count(key)) {
                    Simulator::Cancel(m_units[appid][key].timer);
                    m_units[appid].erase(key);
                }
            }
            NS_LOG_INFO("Aggregation complete for app " << appid << " key " << key);
            return true;
        }
    }
    return false;
}
void 
UdpAggregator::LogArrivalTime(uint16_t appId, uint32_t key) {
    if (!m_timeLogFile.is_open()) return;
    
    // 只记录 appId 0,1,2 且 key=10 的数据
    if (appId <= 2 && key == 10) {
        double timestamp = Simulator::Now().GetSeconds() * 1000.0;
        timestamp -= 2000.0;  // 减去2000毫秒
        m_timeLogFile << appId << " " << key << " " << timestamp << "\n";
        m_timeLogFile.flush(); // 确保数据立即写入文件
        
        m_cachedSamples[{appId, key}].push_back(timestamp);
        
        std::cout << "Recorded gradient for appId=" << appId 
                  << " key=" << key 
                  << " at time=" << timestamp 
                  << std::endl;
    }
}
void 
UdpAggregator::LoadCachedSamples() {
    std::ifstream file("samples.txt");
    if (!file.is_open()) {
        NS_LOG_ERROR("Failed to open samples.txt");
        return;
    }

    // Clear existing samples
    m_cachedSamples.clear();

    uint16_t appId;
    uint32_t key;
    double timestamp;
    
    // 读取数据并按应用ID分组
    std::map<uint16_t, std::vector<double>> arrivalTimes;
    
    while (file >> appId >> key >> timestamp) {
        if (appId <= 2 && key == 10) {  // 只处理 appId 0,1,2 且 key=10 的数据
            arrivalTimes[appId].push_back(timestamp);
        }
    }
    file.close();

    // 对每个应用计算lambda和最优T
    for (auto& [appId, times] : arrivalTimes) {
        if (times.size() < 2) continue;

        // 按时间排序
        std::sort(times.begin(), times.end());

        // 计算相邻到达时间间隔
        std::vector<double> intervals;
        for (size_t i = 1; i < times.size(); ++i) {
            intervals.push_back(times[i] - times[i-1]);
        }

        // 计算平均到达间隔（毫秒）
        double avgInterval = std::accumulate(intervals.begin(), intervals.end(), 0.0) / intervals.size();
        
        // 计算到达率 lambda（每毫秒的数据包数）
        double lambda = 1.0 / avgInterval;

        // 计算最优时间窗口 T（毫秒）
        // 使用公式：T = sqrt(2C/λh)，这里假设 C=1（聚合成本）和 h=1（持有成本）
        double optimalT = std::sqrt(2.0 / lambda);

        // 将最优T从毫秒转换为秒
        double optimalT_seconds = optimalT / 1000.0;

        // 更新应用的时间窗口
        m_currentT[appId] = optimalT_seconds;

        NS_LOG_INFO("App " << appId << ": lambda=" << lambda << " packets/ms, optimal T=" 
                   << optimalT << "ms (" << optimalT_seconds << "s)");
        
        std::cout << "App " << appId << ": lambda=" << lambda << " packets/ms, optimal T=" 
                 << optimalT << "ms (" << optimalT_seconds << "s)" << std::endl;
    }
}
void 
UdpAggregator::ForceFlush(uint16_t appid, uint32_t key) {
    auto& appMap = m_units[appid];
    if (appMap.find(key) == appMap.end()) return;

    AggregatorUnit& unit = appMap[key];
    if (unit.isFlushing) return;

    uint32_t index = app_and_key_to_bitmap_index[appid][key];
    if (appid == 1) {
        NS_LOG_WARN("Force flushing app1 key " << key << " with " << count_pkt[index] << "/" << n << " packets");
    }

    upload(unit.mergedPacket, appid, key, 0, 0);
    unit.isFlushing = true;
    app_and_key_forwarded[appid][key] = 1;
    cleanaggregator(appid, key);
}

void 
UdpAggregator::cleanaggregator(uint16_t appid, uint32_t key)
{
    if(!isexist(appid,key)){
        return;
    }
    
    uint32_t index = app_and_key_to_bitmap_index[appid][key];
    
    // 确保所有定时器都被取消
    if (onofftimewindow) {
        m_merged += 1;
        m_merged_packet += count_pkt[index];
        
        if (m_units.count(appid)) {
            auto& app_units = m_units[appid];
            if (app_units.count(key)) {
                Simulator::Cancel(app_units[key].timer);
                app_units.erase(key);
            }
        }
    }

    // 重置所有状态
    count_pkt[index] = 0;
    aggr_collision[index] = 0;
    idx_appid[index] = 65534;
    
    for(uint32_t i=0; i < max_host_num; i++){
        bitmap[index][i] = false;
        interval[index][i] = -1;
    }
    
    app_and_key_to_bitmap_index[appid].erase(key);
    unused.push_back(index);
    start_time_index[index] = -1;

    if (appid == 1) {
        NS_LOG_INFO("Cleaned aggregator for app1 key " << key);
    }
    
    if(key%1300==0){
        NS_LOG_LYJ("APP"<< appid << " RStime "<<Simulator::Now ().As (Time::S) <<" RestAggregator "<<unused.size());
    }
}

void 
UdpAggregator::broadcast(Ptr<Socket> socket,Ptr<Packet> packet,uint16_t appid, uint32_t key, uint8_t is_bpAggr, uint16_t is_Col, uint32_t syn_urgent, uint8_t set_ecn)
{
    SwitchHeader tmp_header;

    tmp_header.SetKey (key); //set the key in header
    //tmp_header.SetHostid(host_id);
    tmp_header.SetACK(1);
    tmp_header.SetAppID(appid);
    tmp_header.SetDelta(syn_urgent);
    // tmp_header.SetDelta(totaldelta);
    // std::cout<<"debug udp-aggregator broadcast collision_marking "<< collision_marking<<std::endl;
    tmp_header.SetCollision(is_Col);
    tmp_header.SetMACK(set_ecn);
    packet->AddHeader(tmp_header);




    // socket->SendTo  (packet, 0, from);
    int numhost = 0;
    for (auto i =group[appid].begin(); i != group[appid].end(); i++){
        SwmlRouteTag swmlroutetag;
        packet->RemovePacketTag(swmlroutetag);
        swmlroutetag.SetSwmlRouteTag(appid, numhost);
        
        packet->AddPacketTag(swmlroutetag);
        numhost++;
        socket->SendTo(packet,0,*i);
          if (InetSocketAddress::IsMatchingType (*i))
            {
              NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " server sent " << packet->GetSize () << " bytes to " <<
                           InetSocketAddress::ConvertFrom (*i).GetIpv4 () << " port " <<
                           InetSocketAddress::ConvertFrom (*i).GetPort ());
              // if(key == 0) std::cout<<"At time " << Simulator::Now ().As (Time::S) << " server sent " << packet->GetSize () << " bytes to " <<
              //              InetSocketAddress::ConvertFrom (*i).GetIpv4 () << " port " <<
              //              InetSocketAddress::ConvertFrom (*i).GetPort ()<<" key "<<key<<std::endl;
            }
    }
}

void 
UdpAggregator::upload(Ptr<Packet> packet,uint16_t appid, uint32_t key, uint8_t is_bpAggr, uint32_t syn_urgent)
{
    SwitchHeader tmp_header;
    tmp_header.SetKey (key); //set the key in header
    //tmp_header.SetHostid(aggregate_id);
    tmp_header.SetACK(0);
    tmp_header.SetCollision(0);
    // tmp_header.SetDelta(totaldelta);
    tmp_header.SetAppID(appid);
    tmp_header.SetbpAggr(is_bpAggr);
    tmp_header.SetDelta(syn_urgent);
    uint32_t index = app_and_key_to_bitmap_index[appid][key];
    tmp_header.SetTotal(count_pkt[index]);
    tmp_header.SetMerged(1);

    packet->AddHeader(tmp_header);
    // socket->SendTo  (packet, 0, from);
    NS_LOG_INFO("Before send upload ["<<appid<<", "<<key<<"]");
    m_socket_for_up[appid]->Send(packet);
    NS_LOG_INFO("After upload ["<<appid<<", "<<key<<"]");
    if (Ipv4Address::IsMatchingType (multiple_address_for_up[appid]))
        {
        NS_LOG_INFO ("Ipv4Address::IsMatchingType At time " << Simulator::Now ().As (Time::S) << " aggergator sent "<< packet->GetSize () <<" to " <<
                    Ipv4Address::ConvertFrom (multiple_address_for_up[appid]) << " port " << multiple_port_for_up[appid]);
        }
    else if (InetSocketAddress::IsMatchingType (multiple_address_for_up[appid]))
        {
        NS_LOG_INFO ("InetSocketAddress::IsMatchingType At time " << Simulator::Now ().As (Time::S) << " aggergator sent  to " <<
                    InetSocketAddress::ConvertFrom (multiple_address_for_up[appid]).GetIpv4 () << " port " << 
                    //InetSocketAddress::ConvertFrom (m_peerAddress).GetPort ());
                    9);
        }
}


//not used
void 
UdpAggregator::HandleReadForUp (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);

  Ptr<Packet> packet;
  Address from;
  Address localAddress;
  while ((packet = socket->RecvFrom (from)))
    {
      socket->GetSockName (localAddress);
      m_rxTrace (packet);
      m_rxTraceWithAddresses (packet, from, localAddress);
  

      SwitchHeader tmp_header;
      packet->RemoveHeader (tmp_header);
      uint32_t recv_key = tmp_header.GetKey();
      //uint32_t recv_id = tmp_header.GetHostid();
      uint16_t isAck = tmp_header.GetACK();
      //uint16_t isCollision = tmp_header.GetCollision();
      //uint16_t isEnd = tmp_header.GetEnd();
      uint16_t app_id = tmp_header.GetAppID();

      if (InetSocketAddress::IsMatchingType (from))
            {
            NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " server received " << packet->GetSize () << " bytes from " <<
                        InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
                        InetSocketAddress::ConvertFrom (from).GetPort ());
            }
      if(isAck==1){
        NS_LOG_INFO("ack");
        //broadcast(socket,packet,app_id,recv_key);
        
      }
    }
}

void 
UdpAggregator::outputthrought()
{

  std::cout<< "AGGR TIME " << Simulator::Now ().As (Time::S) << " minutil " << 1 - 1.0 * minspace / pool_size << std::endl;
  // for(uint32_t i = 0; i < max_host_num; i++){
  //       double rates =  (double )m_count_sent[i] * 256 /1000/1000* 8*1000;
  //       if(rates >0)
  //         std::cout<<"host "<<i<<","<< rates  <<"Mbps"<<std::endl;
  //        m_count_sent[i] = 0;
         
  //   }
  std::cout<< "AGGR TIME " << Simulator::Now ().As (Time::S) << " minutil " << 1 - 1.0 * minspace / pool_size << std::endl;
  //  std::cout << "total= " << m_merged+m_forwarded
  //   << " TotalMerged=" << m_merged_packet
  //   << " Merged=" << m_merged
  //   << " Forwarded=" << m_forwarded 
  //   << " TotolMerged-Merged=" << m_merged_packet - m_merged
  //   << " timeoutForward=" << m_timeout_forwarded
  //   << "\n";
    // std::cout << "start_time " << start_time.As(Time::NS)  <<  std::endl;
    std::vector<int> cnt(10);
    int empty = 0;
    for (int i = 0; i < pool_size; i ++) {
      if (idx_appid[i] < 10)
        cnt[idx_appid[i]] += 1;
      else  
        empty += 1;
    }
    std::cout << "empty:" << empty << "\napp0:" << cnt[0] << "\napp1:" << cnt[1] << "\napp2:" << cnt[2] << std::endl;
   Simulator::Schedule (Seconds(1), &UdpAggregator::outputthrought, this);
 
  minspace = 0;
  // Simulator::Schedule (Seconds(0.001), &UdpAggregator::outputthrought, this);
}

bool 
UdpAggregator::hasForwarded(uint16_t appid, uint32_t key){
 auto it  = app_and_key_forwarded.find(appid);
  if(it == app_and_key_forwarded.end()){ //not found
    return false;
  }else{
    auto keys_map = it->second;
    auto it2 = keys_map.find(key);
    if(it2 == keys_map.end()){ //not found
      return false;
    }
    else{
      return true; //find it!
    }
  }
}

bool 
UdpAggregator::iskicked(uint16_t appid, uint32_t key){
  auto it  = app_and_key_to_kicked.find(appid);
  if(it == app_and_key_to_kicked.end()){ //not found
    return false;
  }else{
    auto keys_map = it->second;
    auto it2 = keys_map.find(key);
    if(it2 == keys_map.end()){ //not found
      return false;
    }
    else{
      return true; //find it!
    }
  }
}

bool 
UdpAggregator::isexist(uint16_t appid, uint32_t key){
  auto it  = app_and_key_to_bitmap_index.find(appid);
  if(it == app_and_key_to_bitmap_index.end()){ //not found
    return false;
  }else{
    auto keys_map = it->second;
    auto it2 = keys_map.find(key);
    if(it2 == keys_map.end()){ //not found
      return false;
    }
    else{
      return true; //find it!
    }
  }
}

bool 
UdpAggregator::updateindexmap(uint16_t appid, uint32_t key, uint32_t index){
  app_and_key_to_bitmap_index[appid][key] = index;
  return true;
}

// void 
// UdpAggregator::setremotes(std::vector<Address> address, std::vector<uint16_t> port){
//     for(auto it=address.begin(); it!=address.end();it++){
//       multiple_address_for_up.push_back(*it);
//     }
//     for(auto it=port.begin(); it!=port.end();it++){
//       multiple_port_for_up.push_back(*it);
//     }
// }

void 
UdpAggregator::setremotes(std::vector<Address> address, std::vector<uint16_t> port, uint32_t appid){
    
    for(auto it=address.begin(); it!=address.end();it++){
      multiple_address_for_up[appid]=(*it);
      std::cout<<"debug udp-aggregator.cc setremotes address "<<multiple_address_for_up[appid]<<" appid "<<appid<<std::endl;
    }
    for(auto it=port.begin(); it!=port.end();it++){
      multiple_port_for_up[appid]=(*it);
      std::cout<<"debug udp-aggregator.cc setremotes port "<<multiple_port_for_up[appid]<<" appid "<<appid<<std::endl;
    }
}


} // Namespace ns3