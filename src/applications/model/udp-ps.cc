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

#include "udp-ps.h"
#include "udp_switchml_header.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("UdpPsApplication");

NS_OBJECT_ENSURE_REGISTERED (UdpPs);

TypeId
UdpPs::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::UdpPs")
    .SetParent<Application> ()
    .SetGroupName("Applications")
    .AddConstructor<UdpPs> ()
    .AddAttribute ("Port", "Port on which we listen for incoming packets.",
                   UintegerValue (9),
                   MakeUintegerAccessor (&UdpPs::m_port),
                   MakeUintegerChecker<uint16_t> ())
    .AddTraceSource ("Rx", "A packet has been received",
                     MakeTraceSourceAccessor (&UdpPs::m_rxTrace),
                     "ns3::Packet::TracedCallback")
    .AddTraceSource ("RxWithAddresses", "A packet has been received",
                     MakeTraceSourceAccessor (&UdpPs::m_rxTraceWithAddresses),
                     "ns3::Packet::TwoAddressTracedCallback")
    .AddAttribute ("Toalworker", "Toal worker",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpPs::total_worker),
                   MakeUintegerChecker<uint32_t> ())
  ;
  return tid;
}

UdpPs::UdpPs ()
{
  NS_LOG_FUNCTION (this);
}

UdpPs::~UdpPs()
{
  NS_LOG_FUNCTION (this);
  m_socket = 0;
  m_socket6 = 0;
}

void
UdpPs::DoDispose (void)
{
  NS_LOG_FUNCTION (this);
  Application::DoDispose ();
}

void 
UdpPs::StartApplication (void)
{
  NS_LOG_FUNCTION (this);

  if (m_socket == 0)
    {
      TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
      m_socket = Socket::CreateSocket (GetNode (), tid);
      InetSocketAddress local = InetSocketAddress (Ipv4Address::GetAny (), m_port);
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

  m_socket->SetRecvCallback (MakeCallback (&UdpPs::HandleRead, this));
  m_socket6->SetRecvCallback (MakeCallback (&UdpPs::HandleRead, this));
  m_count_sent=0;
  m_received = 0;
  m_totalMerged = 0;
  mergetimes = 0;
  mx = 0;
  //outputthrought();
  Simulator::Schedule(Seconds(1.0), &UdpPs::outputthrought, this);
}
void 
UdpPs::StopApplication ()
{
  NS_LOG_FUNCTION (this);

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
UdpPs::HandleRead (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);
  //NS_LOG_LYJ("rece in ps!");
  Ptr<Packet> packet;
  Address from;
  Address localAddress;
  while ((packet = socket->RecvFrom (from)))
    {
      // std::cout<<"debug udpps recv packet"<<std::endl;
      socket->GetSockName (localAddress);
      m_rxTrace (packet);
      m_rxTraceWithAddresses (packet, from, localAddress);
      

      // packet->RemoveAllPacketTags ();
      packet->RemoveAllByteTags ();

      NS_LOG_INFO ("Echoing packet in ps");
      SwitchHeader tmp_header;
      packet->RemoveHeader (tmp_header);
      uint16_t appid = tmp_header.GetAppID();
      uint32_t recv_key = tmp_header.GetKey();
      uint16_t isCollision = tmp_header.GetCollision();
      uint8_t is_bpAggr = tmp_header.GetbpAggr();
      uint32_t syn_urgent = tmp_header.GetDelta();
      uint16_t total = tmp_header.GetTotal();
      uint8_t isMerged = tmp_header.GetMerged();
      m_totalReceived += total;
      if (isMerged) {
        merge[appid][recv_key] = 1;
        m_totalMerged += total;
        mergetimes += 1;
        m_received += 1;
      } else {

      mp[appid][recv_key] += 1;
        m_received += total;
      }
      if (m_totalReceived > mx) {
        mx = m_totalReceived;
        totalTime = Simulator::Now();
      }
      /*
      UDPecn tag;
      if(packet->RemovePacketTag (tag)){
        NS_LOG_LYJ("get direct pkt with ecn: "<<(uint32_t)tag.GetECN());
      }
      */
      // if (appid == 1) printf("PS debug iscol %d isbp %d seq %d\n", isCollision, is_bpAggr, recv_key);
      if (isCollision==1 || is_bpAggr){
        
        NS_LOG_INFO("get udpps collision pkt key is "<<recv_key);
        UDPecn tag;
        if(packet->RemovePacketTag (tag)){
          NS_LOG_INFO("get direct pkt with ecn: "<<(uint32_t)tag.GetECN());
        }
        update_ecn(recv_key,tag.GetECN());

        // std::cout<<"debug recv iscollision"<<std::endl;

        if(aggregate_pkt(recv_key)){
          NS_LOG_INFO("aggregate finish key is "<<recv_key);

          m_count_sent++;
          tmp_header.SetACK(1);
   
          tmp_header.SetCollision(isCollision);
          tmp_header.SetDelta(syn_urgent);
          tmp_header.SetbpAggr(is_bpAggr);
          tmp_header.SetMACK(key_to_ecn_mark[recv_key]);//for a2tp
          packet->AddHeader(tmp_header);

          // for cc
          // lyj add
          tag.SetECN (key_to_ecn_mark[recv_key]);
          if(key_to_ecn_mark[recv_key]){
            // std::cout<<"debug udp-ps collision recv_ecn key "<<recv_key<<std::endl;
            //NS_LOG_LYJ("reduce the rate");
          }
          key_to_ecn_mark.erase(recv_key);
          // tag.SetECN(0); //debug
          packet->AddPacketTag (tag);
          socket->SendTo (packet, 0, from);
        }
      }
      else{
        // std::cout<<"debug udpps recv aggr"<<std::endl;
        tmp_header.SetACK(1);

        tmp_header.SetCollision(0);
        tmp_header.SetDelta(syn_urgent);
        UDPecn tag;
        if(packet->RemovePacketTag (tag)){
          NS_LOG_INFO("get direct pkt with ecn: "<<(uint32_t)tag.GetECN());
        }
        if (tag.GetECN()){
          // std::cout<<"debug udp-ps aggr recv_ecn key "<<recv_key<<std::endl;
          ;
        }
        tmp_header.SetMACK(tag.GetECN()); //for a2tp
        packet->AddPacketTag (tag);
        m_count_sent++;
        //NS_LOG_LYJ("aggregate finish key is "<<recv_key);
        packet->AddHeader(tmp_header);
        socket->SendTo (packet, 0, from);
        if (InetSocketAddress::IsMatchingType (from))
        {
          NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " server sent " << packet->GetSize () << " bytes to " <<
                       InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
                       InetSocketAddress::ConvertFrom (from).GetPort ());
        }
      }
      
      if (Inet6SocketAddress::IsMatchingType (from))
        {
          NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " server sent " << packet->GetSize () << " bytes to " <<
                       Inet6SocketAddress::ConvertFrom (from).GetIpv6 () << " port " <<
                       Inet6SocketAddress::ConvertFrom (from).GetPort ());
        }
    }
}

bool
UdpPs::aggregate_pkt(uint32_t key)
{
    if(key_to_bitmap_index.find(key)==key_to_bitmap_index.end()){
        key_to_bitmap_index.insert(std::pair<uint32_t, uint32_t>(key,0));
        //get time
        //start_time_index[index] = Simulator::Now ().GetMicroSeconds(); 
    }
    key_to_bitmap_index.find(key)->second = key_to_bitmap_index.find(key)->second+1;
    if(key_to_bitmap_index.find(key)->second==total_worker){
      key_to_bitmap_index.erase(key);
      return true;
    }
    else{
      return false;
    }
}

void
UdpPs::update_ecn(uint32_t key, uint8_t ecnmark)
{
    if(key_to_ecn_mark.find(key)==key_to_ecn_mark.end()){
        key_to_ecn_mark.insert(std::pair<uint32_t, uint32_t>(key, 0));
        //get time
        //start_time_index[index] = Simulator::Now ().GetMicroSeconds(); 
    }
    if(ecnmark!=0){
      key_to_ecn_mark[key] = (uint32_t)ecnmark;
    }
}

void 
UdpPs::outputthrought()
{
  std::cout << "[" << Simulator::Now().As(Time::S) 
  << "] Total Received: " << m_received
  << " Merged Sum: " << m_totalMerged 
  << " mergetimes:" << mergetimes 
  << " ans:" << m_totalMerged - mergetimes 
  << std::endl;

  std::cout << "finish_time " << totalTime.As(Time::S) << " total="<< mx <<  std::endl;
  // std::cout << mp[0].size() << " " << mp[1].size() << " " << mp[2].size() << '\n';

  // std::cout << merge[0].size() << " " << merge[1].size() << " " << merge[2].size() << '\n';
  // for (auto x : mp[0]) if (x.second != 8) std:: cout << x.second << " ";
  Simulator::Schedule(Seconds(1.0), &UdpPs::outputthrought, this);
  // double rates =  (double )m_count_sent * 256 /1000/1000* 8*10000;
  // //std::cout<<"count"<<m_count_sent<<std::endl;
  // if(rates>0) std::cout<<"thrtime "<<Simulator::Now ().As (Time::S)<< " aggreate throughput "<< rates  <<" Mbps "<<std::endl;
  // m_count_sent = 0;
  //Simulator::Schedule (Seconds(0.0001), &UdpPs::outputthrought, this);
}

} // Namespace ns3
