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
#include "ns3/nstime.h"
#include "ns3/inet-socket-address.h"
#include "ns3/inet6-socket-address.h"
#include "ns3/socket.h"
#include "ns3/simulator.h"
#include "ns3/socket-factory.h"
#include "ns3/packet.h"
#include "ns3/uinteger.h"
#include "ns3/trace-source-accessor.h"
#include "udp-worker.h"
#include "udp_switchml_header.h"
#include "ns3/double.h"
#include <math.h>


namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("UdpWorkerApplication");

NS_OBJECT_ENSURE_REGISTERED (UdpWorker);

TypeId
UdpWorker::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::UdpWorker")
    .SetParent<Application> ()
    .SetGroupName("Applications")
    .AddConstructor<UdpWorker> ()
    .AddAttribute ("MaxPackets", 
                   "The maximum number of packets the application will send",
                   UintegerValue (100),
                   MakeUintegerAccessor (&UdpWorker::m_count),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("Interval", 
                   "The time to wait between packets",
                   TimeValue (Seconds (1.0)),
                   MakeTimeAccessor (&UdpWorker::m_interval),
                   MakeTimeChecker ())
    .AddAttribute ("RemoteAddress", 
                   "The destination Address of the outbound packets",
                   AddressValue (),
                   MakeAddressAccessor (&UdpWorker::m_peerAddress),
                   MakeAddressChecker ())
    .AddAttribute ("RemotePort", 
                   "The destination port of the outbound packets",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpWorker::m_peerPort),
                   MakeUintegerChecker<uint16_t> ())
    .AddAttribute ("PacketSize", "Size of echo data in outbound packets",
                   UintegerValue (100),
                   MakeUintegerAccessor (&UdpWorker::SetDataSize,
                                         &UdpWorker::GetDataSize),
                   MakeUintegerChecker<uint32_t> ())
    .AddTraceSource ("Tx", "A new packet is created and is sent",
                     MakeTraceSourceAccessor (&UdpWorker::m_txTrace),
                     "ns3::Packet::TracedCallback")
    .AddTraceSource ("Rx", "A packet has been received",
                     MakeTraceSourceAccessor (&UdpWorker::m_rxTrace),
                     "ns3::Packet::TracedCallback")
    .AddTraceSource ("TxWithAddresses", "A new packet is created and is sent",
                     MakeTraceSourceAccessor (&UdpWorker::m_txTraceWithAddresses),
                     "ns3::Packet::TwoAddressTracedCallback")
    .AddTraceSource ("RxWithAddresses", "A packet has been received",
                     MakeTraceSourceAccessor (&UdpWorker::m_rxTraceWithAddresses),
                     "ns3::Packet::TwoAddressTracedCallback")
    .AddAttribute ("Maxbytes", 
                   "The maximum bytes of packets the application will send",
                   UintegerValue (100000),
                   MakeUintegerAccessor (&UdpWorker::m_total_size),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("CWND", 
                   "The send windows of the application will send",
                   UintegerValue (100),
                   MakeUintegerAccessor (&UdpWorker::m_cwnd),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("MAXCWND", 
                   "The send max windows of the application will send",
                   UintegerValue (650),
                   MakeUintegerAccessor (&UdpWorker::m_cwndmax), 
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("Host_number", "host number",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpWorker::m_host_num),
                   MakeUintegerChecker<uint8_t> ())
    .AddAttribute ("HOSTID", 
                   "unique host id",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpWorker::m_host_id),
                   MakeUintegerChecker<uint8_t> ())
    .AddAttribute ("APPID", 
                   "unique app id",
                   UintegerValue (0),
                   MakeUintegerAccessor (&UdpWorker::m_appid),
                   MakeUintegerChecker<uint16_t> ())
    .AddAttribute ("ONOFFCC", 
                   "cc enbale",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpWorker::onoffcc),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("ONOFFASYCC", 
                   "asycc enbale",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpWorker::onoffasycc),
                   MakeUintegerChecker<uint32_t> ())
    .AddAttribute ("ONOFFLZYCC", 
                   "asycc enbale",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpWorker::onofflzycc),
                   MakeUintegerChecker<uint32_t> ())

    .AddAttribute ("ONOFFAWNDCC", 
                   "awndcc enbale",
                   UintegerValue (1),
                   MakeUintegerAccessor (&UdpWorker::onoffawndcc),
                   MakeUintegerChecker<uint32_t> ())                   

    .AddAttribute ("USEDAGGR", 
                   "the size of used aggregators",
                   UintegerValue (40000),
                   MakeUintegerAccessor (&UdpWorker::used_aggr),
                   MakeUintegerChecker<uint32_t> ())      
      .AddAttribute ("BASERTT", 
                   "basertt",
                   DoubleValue (40000),
                   MakeDoubleAccessor (&UdpWorker::baseRTT),
                   MakeDoubleChecker<double> ())               
  ;
  return tid;
}

UdpWorker::UdpWorker ()
{
  NS_LOG_FUNCTION (this);
  m_sent = 0;
  m_socket = 0;
  m_sendEvent = EventId ();
  m_data = 0;
  m_dataSize = 0;
  m_cc = simplecc();
  curseq_ = 0;
  max_ack = 0;
  last_ack = 0;
  ack_bitmap.clear();
  last_awnd_time = 0;
  last_urgent_time = 0;
  updateurgent = 1;
  baseRTT = 0.000016;
}

UdpWorker::~UdpWorker()
{
  NS_LOG_FUNCTION (this);
  m_socket = 0;

  delete [] m_data;
  m_data = 0;
  m_dataSize = 0;
}

void 
UdpWorker::SetRemote (Address ip, uint16_t port)
{
  NS_LOG_FUNCTION (this << ip << port);
  m_peerAddress = ip;
  m_peerPort = port;
}

void 
UdpWorker::SetRemote (Address addr)
{
  NS_LOG_FUNCTION (this << addr);
  m_peerAddress = addr;
}

void
UdpWorker::DoDispose (void)
{
  NS_LOG_FUNCTION (this);
  Application::DoDispose ();
}

void 
UdpWorker::StartApplication (void)
{
  NS_LOG_FUNCTION (this);

  if (m_socket == 0)
    {
      TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
      m_socket = Socket::CreateSocket (GetNode (), tid);
      if (Ipv4Address::IsMatchingType(m_peerAddress) == true)
        {
          if (m_socket->Bind () == -1)
            {
              NS_FATAL_ERROR ("Failed to bind socket");
            }
          m_socket->Connect (InetSocketAddress (Ipv4Address::ConvertFrom(m_peerAddress), m_peerPort));
        }
      else if (Ipv6Address::IsMatchingType(m_peerAddress) == true)
        {
          if (m_socket->Bind6 () == -1)
            {
              NS_FATAL_ERROR ("Failed to bind socket");
            }
          m_socket->Connect (Inet6SocketAddress (Ipv6Address::ConvertFrom(m_peerAddress), m_peerPort));
        }
      else if (InetSocketAddress::IsMatchingType (m_peerAddress) == true)
        {
          if (m_socket->Bind () == -1)
            {
              NS_FATAL_ERROR ("Failed to bind socket");
            }
          m_socket->Connect (m_peerAddress);
        }
      else if (Inet6SocketAddress::IsMatchingType (m_peerAddress) == true)
        {
          if (m_socket->Bind6 () == -1)
            {
              NS_FATAL_ERROR ("Failed to bind socket");
            }
          m_socket->Connect (m_peerAddress);
        }
      else
        {
          NS_ASSERT_MSG (false, "Incompatible address type: " << m_peerAddress);
        }
    }

  m_socket->SetRecvCallback (MakeCallback (&UdpWorker::SwitchmlHandleRead, this));
  m_socket->SetAllowBroadcast (true);
  setindex();
  ScheduleTransmit (Seconds (0.));
  max_ack = 0;
  last_ack=0;
  last_mack = 0;
  ack_bitmap.clear();
  logthroughput(); //for large scale simulation

  Simulator::Schedule(Seconds(1.0), &UdpWorker::outputthrought, this);
}

void 
UdpWorker::StopApplication ()
{
  NS_LOG_FUNCTION (this);

  if (m_socket != 0) 
    {
      m_socket->Close ();
      m_socket->SetRecvCallback (MakeNullCallback<void, Ptr<Socket> > ());
      m_socket = 0;
    }

  Simulator::Cancel (m_sendEvent);
}

void 
UdpWorker::SetDataSize (uint32_t dataSize)
{
  NS_LOG_FUNCTION (this << dataSize);

  //
  // If the client is setting the echo packet data size this way, we infer
  // that she doesn't care about the contents of the packet at all, so 
  // neither will we.
  //
  delete [] m_data;
  m_data = 0;
  m_dataSize = 0;
  m_size = dataSize;
}

uint32_t 
UdpWorker::GetDataSize (void) const
{
  NS_LOG_FUNCTION (this);
  return m_size;
}

void 
UdpWorker::SetFill (std::string fill)
{
  NS_LOG_FUNCTION (this << fill);

  uint32_t dataSize = fill.size () + 1;

  if (dataSize != m_dataSize)
    {
      delete [] m_data;
      m_data = new uint8_t [dataSize];
      m_dataSize = dataSize;
    }

  memcpy (m_data, fill.c_str (), dataSize);

  //
  // Overwrite packet size attribute.
  //
  m_size = dataSize;
}

void 
UdpWorker::SetFill (uint8_t fill, uint32_t dataSize)
{
  NS_LOG_FUNCTION (this << fill << dataSize);
  if (dataSize != m_dataSize)
    {
      delete [] m_data;
      m_data = new uint8_t [dataSize];
      m_dataSize = dataSize;
    }

  memset (m_data, fill, dataSize);

  //
  // Overwrite packet size attribute.
  //
  m_size = dataSize;
}

void 
UdpWorker::SetFill (uint8_t *fill, uint32_t fillSize, uint32_t dataSize)
{
  NS_LOG_FUNCTION (this << fill << fillSize << dataSize);
  if (dataSize != m_dataSize)
    {
      delete [] m_data;
      m_data = new uint8_t [dataSize];
      m_dataSize = dataSize;
    }

  if (fillSize >= dataSize)
    {
      memcpy (m_data, fill, dataSize);
      m_size = dataSize;
      return;
    }

  //
  // Do all but the final fill.
  //
  uint32_t filled = 0;
  while (filled + fillSize < dataSize)
    {
      memcpy (&m_data[filled], fill, fillSize);
      filled += fillSize;
    }

  //
  // Last fill may be partial
  //
  memcpy (&m_data[filled], fill, dataSize - filled);

  //
  // Overwrite packet size attribute.
  //
  m_size = dataSize;
}

void 
UdpWorker::ScheduleTransmit (Time dt)
{
  NS_LOG_FUNCTION (this << dt);
  m_sendEvent = Simulator::Schedule (dt, &UdpWorker::SwitchmlSend, this);
  //Simulator::Schedule (Seconds(1), &UdpWorker::logthroughput, this);
}

void 
UdpWorker::Send (void)
{
  NS_LOG_FUNCTION (this);

  NS_ASSERT (m_sendEvent.IsExpired ());

  Ptr<Packet> p;
  if (m_dataSize)
    {
      //
      // If m_dataSize is non-zero, we have a data buffer of the same size that we
      // are expected to copy and send.  This state of affairs is created if one of
      // the Fill functions is called.  In this case, m_size must have been set
      // to agree with m_dataSize
      //
      NS_ASSERT_MSG (m_dataSize == m_size, "UdpWorker::Send(): m_size and m_dataSize inconsistent");
      NS_ASSERT_MSG (m_data, "UdpWorker::Send(): m_dataSize but no m_data");
      p = Create<Packet> (m_data, m_dataSize);
    }
  else
    {
      //
      // If m_dataSize is zero, the client has indicated that it doesn't care
      // about the data itself either by specifying the data size by setting
      // the corresponding attribute or by not calling a SetFill function.  In
      // this case, we don't worry about it either.  But we do allow m_size
      // to have a value different from the (zero) m_dataSize.
      //
      p = Create<Packet> (m_size);
    }
  
  SwitchHeader tmp_header;
  
  UDPecn tag;
  tag.SetECN (0);
  p->AddPacketTag (tag);

  SwmlRouteTag routetag;
  routetag.SetSwmlRouteTag(m_appid, m_host_id);
  p->AddPacketTag(routetag);

  tmp_header.SetKey (2);
  p->AddHeader(tmp_header);
  Address localAddress;
  m_socket->GetSockName (localAddress);
  // call to the trace sinks before the packet is actually sent,
  // so that tags added to the packet can be sent as well
  m_txTrace (p);
  if (Ipv4Address::IsMatchingType (m_peerAddress))
    {
      m_txTraceWithAddresses (p, localAddress, InetSocketAddress (Ipv4Address::ConvertFrom (m_peerAddress), m_peerPort));
    }
  else if (Ipv6Address::IsMatchingType (m_peerAddress))
    {
      m_txTraceWithAddresses (p, localAddress, Inet6SocketAddress (Ipv6Address::ConvertFrom (m_peerAddress), m_peerPort));
    }
  m_socket->Send (p);
  ++m_sent;

  if (Ipv4Address::IsMatchingType (m_peerAddress))
    {
      NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                   Ipv4Address::ConvertFrom (m_peerAddress) << " port " << m_peerPort);
    }
  else if (Ipv6Address::IsMatchingType (m_peerAddress))
    {
      NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                   Ipv6Address::ConvertFrom (m_peerAddress) << " port " << m_peerPort);
    }
  else if (InetSocketAddress::IsMatchingType (m_peerAddress))
    {
      NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                   InetSocketAddress::ConvertFrom (m_peerAddress).GetIpv4 () << " port " << InetSocketAddress::ConvertFrom (m_peerAddress).GetPort ());
    }
  else if (Inet6SocketAddress::IsMatchingType (m_peerAddress))
    {
      NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                   Inet6SocketAddress::ConvertFrom (m_peerAddress).GetIpv6 () << " port " << Inet6SocketAddress::ConvertFrom (m_peerAddress).GetPort ());
    }

  if (m_sent < m_count) 
    {
      ScheduleTransmit (m_interval);
    }
}

void
UdpWorker::HandleRead (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);
  Ptr<Packet> packet;
  Address from;
  Address localAddress;
  while ((packet = socket->RecvFrom (from)))
    {
      if (InetSocketAddress::IsMatchingType (from))
        {
          NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client received " << packet->GetSize () << " bytes from " <<
                       InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
                       InetSocketAddress::ConvertFrom (from).GetPort ());
          // std::cout << "lyj:At time " << Simulator::Now ().As (Time::S) << " client received " << packet->GetSize () << " bytes from " <<
          //              InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
          //              InetSocketAddress::ConvertFrom (from).GetPort ()<<std::endl;
        }
      else if (Inet6SocketAddress::IsMatchingType (from))
        {
          NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client received " << packet->GetSize () << " bytes from " <<
                       Inet6SocketAddress::ConvertFrom (from).GetIpv6 () << " port " <<
                       Inet6SocketAddress::ConvertFrom (from).GetPort ());
        }
      socket->GetSockName (localAddress);
      m_rxTrace (packet);
      m_rxTraceWithAddresses (packet, from, localAddress);
    }
}

//lyj add
void 
UdpWorker::SendPacket (uint32_t index, uint8_t host_id, uint32_t bp_Aggr = 0, uint32_t syn_urgent = 100000, uint8_t set_ecn = 0)
{
  Ptr<Packet> p;
    if (m_dataSize)
      {
        //
        // If m_dataSize is non-zero, we have a data buffer of the same size that we
        // are expected to copy and send.  This state of affairs is created if one of
        // the Fill functions is called.  In this case, m_size must have been set
        // to agree with m_dataSize
        //
        NS_ASSERT_MSG (m_dataSize == m_size, "UdpWorker::Send(): m_size and m_dataSize inconsistent");
        NS_ASSERT_MSG (m_data, "UdpWorker::Send(): m_dataSize but no m_data");
        p = Create<Packet> (m_data, m_dataSize);
      }
    else
      {
        //
        // If m_dataSize is zero, the client has indicated that it doesn't care
        // about the data itself either by specifying the data size by setting
        // the corresponding attribute or by not calling a SetFill function.  In
        // this case, we don't worry about it either.  But we do allow m_size
        // to have a value different from the (zero) m_dataSize.
        //
        p = Create<Packet> (m_size);
      }
    
    SwitchHeader tmp_header;

    UDPecn tag;
    tag.SetECN (0);
    p->AddPacketTag (tag);

    SwmlRouteTag routetag;
    routetag.SetSwmlRouteTag(m_appid, m_host_id);
    p->AddPacketTag(routetag);


    tmp_header.SetKey (index); //set the key in header
    tmp_header.SetHostid(host_id);
    if(host_id!=m_host_id){
      NS_LOG_INFO("i help send "<<index);
    }
    tmp_header.SetACK(0);
    tmp_header.SetCollision(0);
    tmp_header.SetMACK(set_ecn); //for a2tp
    tmp_header.SetAppID(m_appid);
    tmp_header.SetHostnum(m_host_num);
    tmp_header.SetbpAggr(bp_Aggr);
    tmp_header.SetDelta(syn_urgent);
    
    // if(index==last_cwnd_pkt){
    //   tmp_header.SetEnd(1);
    // }
    p->AddHeader(tmp_header);

    //NS_LOG_LYJ("set: "<<(uint32_t)ipTosTag.GetECN());
    Address localAddress;
    m_socket->GetSockName (localAddress);
    // call to the trace sinks before the packet is actually sent,
    // so that tags added to the packet can be sent as well
    m_txTrace (p);
    if (Ipv4Address::IsMatchingType (m_peerAddress))
      {
        m_txTraceWithAddresses (p, localAddress, InetSocketAddress (Ipv4Address::ConvertFrom (m_peerAddress), m_peerPort));
      }
    else if (Ipv6Address::IsMatchingType (m_peerAddress))
      {
        m_txTraceWithAddresses (p, localAddress, Inet6SocketAddress (Ipv6Address::ConvertFrom (m_peerAddress), m_peerPort));
      }
    uint32_t send_data = m_socket->Send (p);
    if(send_data < m_size){
      NS_LOG_INFO("overflow!");
    }
    ++m_sent;

    if (Ipv4Address::IsMatchingType (m_peerAddress))
      {
        NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                    Ipv4Address::ConvertFrom (m_peerAddress) << " port " << m_peerPort);
      }
    else if (Ipv6Address::IsMatchingType (m_peerAddress))
      {
        NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                    Ipv6Address::ConvertFrom (m_peerAddress) << " port " << m_peerPort);
      }
    else if (InetSocketAddress::IsMatchingType (m_peerAddress))
      {
        NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                    InetSocketAddress::ConvertFrom (m_peerAddress).GetIpv4 () << " port " << InetSocketAddress::ConvertFrom (m_peerAddress).GetPort ());
      }
    else if (Inet6SocketAddress::IsMatchingType (m_peerAddress))
      {
        NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client sent " << m_size << " bytes to " <<
                    Inet6SocketAddress::ConvertFrom (m_peerAddress).GetIpv6 () << " port " << Inet6SocketAddress::ConvertFrom (m_peerAddress).GetPort ());
      }
}


uint32_t UdpWorker::UpdateNewLastACK(uint32_t ack_idx){
  ack_bitmap[ack_idx] = true;
  uint32_t new_last_ack = 0;
  uint32_t cur_idx = last_ack;
  // std::cout<<"debug UpdateNewLastACK ack_idx "<<ack_idx<<" cur_idx "<<cur_idx<<" last_ack "<<last_ack<<" max_ack "<<max_ack<<std::endl;
  while (ack_bitmap.find(cur_idx) != ack_bitmap.end() && cur_idx <= max_ack){
    new_last_ack = cur_idx;
    cur_idx ++;
  }
  last_ack_time = Simulator::Now();
  return new_last_ack;
}

void 
UdpWorker::SwitchmlSend (void)
{
  NS_LOG_FUNCTION (this);

  NS_ASSERT (m_sendEvent.IsExpired ());
  m_cc.initcwnd(m_cwnd, m_cwndmax, onoffcc, onoffasycc, used_aggr);

  if(m_cwnd < m_max_index){
    last_cwnd_pkt = m_cwnd-1;
  }else{
    last_cwnd_pkt = m_max_index-1;
  }


  for (uint32_t i = 0; i < m_cwnd && i < m_max_index; i++)
  {
    SendPacket(i,m_host_id, 0, 1);
    curseq_ = i + 1;
    NS_LOG_INFO("send key "<<i);
    m_cc.setlast_send(i);
  }
  // std::cout<<"debug SwitchmlSend"<< std::endl;
}

void
UdpWorker::SwitchmlHandleRead (Ptr<Socket> socket)
{
  NS_LOG_FUNCTION (this << socket);
  Ptr<Packet> packet;
  Address from;
  Address localAddress;
  while ((packet = socket->RecvFrom (from)))
    {
      if (InetSocketAddress::IsMatchingType (from))
        {
          // std::cout<<"debug 1 SwitchmlHandleRead"<<std::endl;
          // for saatp
          if(0&&packet->GetSize () !=270){ //is a notify packet
            uint32_t loc = packet->GetSize ()-1;
            NS_LOG_INFO("i receive the notify, last_ack_is "<<last_ack<<" end pkt is "<< last_cwnd_pkt);
            for(uint32_t i = last_ack; i<last_cwnd_pkt+1;i++){
              SendPacket(i,loc);
            }
            NS_LOG_INFO("i send as "<<loc);
          }
          else{
              // lyj add
            SwitchHeader tmp_header;
            packet->RemoveHeader (tmp_header);
             // for cc
            // lyj add
            UDPecn ipTosTag;
            packet->RemovePacketTag (ipTosTag);
            uint8_t link_ecn = ipTosTag.GetECN();
            // if(packet->RemovePacketTag (ipTosTag)){
            //   //NS_LOG_LYJ("i successfully recevied a ecn: "<<(uint32_t)ipTosTag.GetECN());
            //   if(ipTosTag.GetECN()){
            //     // printf("app %d recv ECN\n", m_appid);
            //     m_cc.ecn_count++;
            //   }
            //   m_cc.update_ecn(ipTosTag.GetECN());
            // }

            uint8_t a2tpecn = tmp_header.GetMACK();
            if(a2tpecn) m_cc.ecn_count++;
            m_cc.update_ecn(a2tpecn);
            uint32_t cur_index = tmp_header.GetKey();
            uint8_t isack = tmp_header.GetACK();
            // uint8_t ismack = tmp_header.GetMACK();
            uint16_t iscol = tmp_header.GetCollision();

            if(iscol){
              // printf("col\n");
            }
            // m_cc.urgent = tmp_header.GetDelta()/100000.0; //syn_urgent

            m_cc.urgent = tmp_header.GetDelta()/100000.0; //syn_urgent


            bool enablea2tp = 1;
            if (onoffawndcc){
              if(cur_index > max_ack){
                max_ack = cur_index;
              }
              if(cur_index >= last_ack){
                // if (cur_index > last_ack + 1) printf("Worker debug appid %d host %d recv %d last_ack %d \n", m_appid, m_host_id, cur_index, last_ack);
                last_ack = UpdateNewLastACK(cur_index);
                // if (m_appid == 1) printf("Worker debug host %d recv %d last_ack %d \n", m_host_id, cur_index, last_ack);
                // last_ack = cur_index; // i think there are no reorder and loss

                if(m_host_id == 0){
                    ;//NS_LOG_LYJ("recieve packet last_ack "<<last_ack);
                }
              }
              if(iscol){
                  // printf("appid %d, hostid %d, cur_index %d\n",m_appid, m_host_id, cur_index);
                m_cc.col_count ++;
              }
              m_cc.recv_ack_num++;
              for(uint32_t i = curseq_; i<= last_ack + m_cc.get_cwnd();i++){
                if(i<m_max_index){
                  uint8_t set_ecn = 0;
                  if (link_ecn && !a2tpecn) set_ecn = 1;
                  if (i <= last_ack + m_cc.get_awnd()){
                    SendPacket(i,m_host_id,0,  (uint32_t) (100000 * m_cc.lastdelta), 0);
                  }else{
                    SendPacket(i,m_host_id,1, (uint32_t) (100000 * m_cc.lastdelta), 0);
                  }
                  curseq_++;
                  //NS_LOG_INFO("send key "<<i);
                  NS_LOG_INFO("host"<<m_host_id<<" send packet curseq_ "<<i);
                  
                  //m_cc.setlast_send(i);
                }
              }
              double tmpnowtime = Simulator::Now().ToDouble(Time::S);
              if (tmpnowtime - last_urgent_time > 1.01*baseRTT && updateurgent == 0){
                updateurgent = 1;
                float tmpurgent = 1.0 * m_cc.recv_ack_num / m_cc.get_cwnd();
                // if (m_cc.ecn_count / m_cc.recv_ack_num > 0.3) tmpurgent = 1;
                // tmpurgent  = tmpurgent > 1 ? 1 : tmpurgent;
                
                float uw = 1.0 / 16;
                m_cc.lastdelta = (1-uw) * m_cc.lastdelta + uw * tmpurgent;
                // if (m_appid == 9) printf("Worker debug host %d last_ack %d tmpurgent %f urgent %f recv %d cwnd %d\n", m_host_id, last_ack, tmpurgent, m_cc.urgent, m_cc.recv_ack_num, m_cc.get_cwnd());
                if (m_host_id == 0){
                  // printf("appid %d, curseq_ %d, col_count %d, recv %d, tmpurgent %f, urgent %f, cwnd %d awnd %d\n", m_appid, curseq_, m_cc.col_count, m_cc.recv_ack_num, tmpurgent, m_cc.urgent, m_cc.get_cwnd(), m_cc.get_awnd());
                  ;// printf("appid %d, hostid %d, curseq_ %d, col_count %d, recv %d, tmpbeta %f, delta %d, awnd %d\n", m_appid, m_host_id, curseq_, m_cc.col_count, m_cc.recv_ack_num, tmpbeta, curseq_-last_ack, m_cc.get_awnd());
                }
                // m_cc.recv_ack_num = 1;
                // last_urgent_time = tmpnowtime;
              }
            }else{

              if(cur_index > max_ack){
                max_ack = cur_index;
              }
              if(cur_index >= last_ack){
                last_ack = UpdateNewLastACK(cur_index);
                // last_ack = cur_index; // i think there are no reorder and loss

                if(m_host_id == 0){
                    ;//NS_LOG_LYJ("recieve packet last_ack "<<last_ack);
                }
              }
              for(uint32_t i = curseq_; i<= last_ack + m_cc.get_cwnd();i++){
                if(i<m_max_index){
                  SendPacket(i,m_host_id);
                  curseq_++;
                  //NS_LOG_INFO("send key "<<i);
                  NS_LOG_INFO("host"<<m_host_id<<" send packet curseq_ "<<i);
                  
                  //m_cc.setlast_send(i);
                }
              }
            }
            
            if(curseq_ >= last_cwnd_pkt){
              
              NS_LOG_LYJ("APP "<<m_appid<<" time "<< Simulator::Now ().As (Time::S)<<" host"<<m_host_id<<" cur cwnd is "<<m_cc.get_cwnd()<<" delta "<< curseq_ - last_ack \
                << " weight "<< m_cc.weight);
              // std::cout<<"APP "<<m_appid<<" time "<< Simulator::Now ().As (Time::S)<<" host"<<(uint32_t)m_host_id<<" cur cwnd is "<<m_cc.get_cwnd()<<" delta "<< curseq_ - last_ack \
              //   << " weight "<< m_cc.weight << " awnd "<<m_cc.get_awnd()<<std::endl;
              if (m_host_id == 0){
                ;
                // printf("appid %d, curseq_ %d, ecn %d cwnd %d\n", m_appid, curseq_, m_cc.ecn_count, m_cc.get_cwnd());
              }
              // if (m_appid == 9) printf("Worker debug host %d last_ack %d urgent %f recv %d cwnd %d awnd %d\n", m_host_id, last_ack, m_cc.urgent, m_cc.recv_ack_num, m_cc.get_cwnd(), m_cc.get_awnd());
              m_cc.adjustcwnd(onofflzycc, onoffcc);
              m_cc.clean_ecn();
              last_cwnd_pkt = curseq_ + m_cc.get_cwnd();

              m_cc.recv_ack_num = 1;
              double tmpnowtime1 = Simulator::Now().ToDouble(Time::S);
              last_urgent_time = tmpnowtime1;
              updateurgent = 0;
              
            }
            

            //m_cc.setlast_ack(cur_index); // i think there are no reorder and loss


           
            /*uint32_t left = m_cc.get_send_left_cound();
            uint32_t right = m_cc.get_send_right_cound();
            NS_LOG_LYJ("send from ["<<left<<", "<<right<<")");
            for(uint32_t i = left; i<right;i++){
              if(i<m_max_index){
                SendPacket(i,m_host_id);
                NS_LOG_INFO("send key "<<i);
                m_cc.setlast_send(i);
              }
            }
            */
            if(m_host_id == 0){
              NS_LOG_INFO("curseq_ last_ack cwnd"<<curseq_<<" "<<last_ack<<" "<<m_cc.get_cwnd());
            }
            // uint32_t next_index = cur_index + m_cwnd;
            // if(next_index<m_max_index){
            //   SendPacket(next_index,m_host_id);
            //   NS_LOG_INFO("send key "<<next_index);
            // }
            NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client received " << packet->GetSize () << " bytes from " <<
                        InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
                        InetSocketAddress::ConvertFrom (from).GetPort ());
            // std::cout << "lyj:At time " << Simulator::Now ().As (Time::S) << " client received " << packet->GetSize () << " bytes from " <<
            //              InetSocketAddress::ConvertFrom (from).GetIpv4 () << " port " <<
            //              InetSocketAddress::ConvertFrom (from).GetPort ()<<std::endl;
            if(isack){
              m_complete++;
              m_count_sent++;
              NS_LOG_INFO("compelet"<<m_complete<<" key is"<<cur_index);
            }
            
            if (m_complete == m_max_index && m_host_id == 0){
              NS_LOG_LYJ("Finish "<< m_appid<<" node "<<GetNode ()->GetId()<<" at "<<Simulator::Now ().As (Time::S));
              std::cout<<"Finish "<< m_appid<<" node "<<GetNode ()->GetId()<<" at "<<Simulator::Now ().As (Time::S)<<std::endl;
            }
          }
          
        }
      else if (Inet6SocketAddress::IsMatchingType (from))
        {
          NS_LOG_INFO ("At time " << Simulator::Now ().As (Time::S) << " client received " << packet->GetSize () << " bytes from " <<
                       Inet6SocketAddress::ConvertFrom (from).GetIpv6 () << " port " <<
                       Inet6SocketAddress::ConvertFrom (from).GetPort ());
        }
      socket->GetSockName (localAddress);
      m_rxTrace (packet);
      m_rxTraceWithAddresses (packet, from, localAddress);
    }
}

void 
UdpWorker::setindex (void)
{
  m_size_of_per_pkt = m_size;
  m_max_index = ( m_total_size + m_size_of_per_pkt -1 ) / m_size_of_per_pkt;
  m_complete = 0;
  m_count_sent = 0;
}

void 
UdpWorker::logthroughput (void)
{
  
  if(m_complete < m_max_index && m_host_id == 0){
  // if(m_complete < m_max_index){
    NS_LOG_LYJ("APP "<<m_appid<<" thrtime "<<Simulator::Now ().As (Time::S)<< " aggregate throughput "<< m_count_sent * m_size_of_per_pkt /1000000000.0 * 8 * 1000  <<" Gbps");
    std::cout<<"APP "<<m_appid<<" thrtime "<<Simulator::Now ().As (Time::S)<< " aggregate throughput "<< m_count_sent * m_size_of_per_pkt /1000000000.0 * 8 * 1000  <<" Gbps "\
    << "[("<<m_complete<<"/"<<m_max_index <<"),"<<1.0*m_complete/m_max_index*100.0<<"%]"<<std::endl;
    if (m_complete <= m_max_index){
      Simulator::Schedule (Seconds(0.001), &UdpWorker::logthroughput, this);
    }
    if (m_complete > mx) {
      mx = m_complete;
      totalTime = Simulator::Now();
    }
  }
  m_count_sent = 0;
}
void 
UdpWorker::outputthrought() {
  std::cout << "[" << Simulator::Now().As(Time::S) 
  << "] Total Send: " << m_complete
  << std::endl;

  std::cout << "finish_send_time " << totalTime.As(Time::S) << " total="<< mx <<  std::endl;
  std::cout << "last_ack_time " << last_ack_time.As(Time::S) << " total="<< mx <<  std::endl;
  Simulator::Schedule(Seconds(1.0), &UdpWorker::outputthrought, this);
}

simplecc::simplecc(){
      cwnd = 0;
      isecn = 0;
      cwndmax_ = 680;
      ssthre_ = cwndmax_;
      lastdelta = 1;
      weight = 1.0;
      ecn_count = 0;
      alpha = 0;
      urgent = 1;
      col_count = 0;
      recv_ack_num = 1;
      beta = 0;
      awnd = 0;

    }

simplecc::~simplecc(){
}

void 
simplecc::initcwnd(uint32_t tmp, uint32_t cwndmax, uint32_t onoffcc, uint32_t onoffasycc, uint32_t used_aggr){
  cwndmax_ = cwndmax;
  ssthre_ = cwndmax_;
  cwnd = cwndmax; //for large scale simulation
  isecn = 0;
  awndmax_ = used_aggr;
  awndmin_ = 0;
  awndssthre_ = used_aggr;

  awnd = used_aggr;
  urgent = 1;
  enablecc =onoffcc;
  enableasycc = onoffasycc;
}


void
simplecc::adjustawnd(float tmpbeta){

; //no use

}


void 
simplecc::adjustcwnd(uint32_t onofflzycc, uint32_t onoffcc){
 

  if(onoffcc){
      float g = 1.0/16;
      float tmpalpha = 1.0 * ecn_count / cwnd;
      tmpalpha = tmpalpha > 1 ? 1 : tmpalpha;
      alpha = g * tmpalpha + (1-g) * alpha;



      if(ecn_count > 0){
          cwnd = cwnd * (1 - alpha/2);
          // cwnd = cwnd * (1 - 1.0/2);
          // cwnd = cwnd * (1 - pow( alpha, 0.25)/2);
          
          ssthre_ = cwnd;
      }else{
        if(cwnd >= ssthre_){
          cwnd = cwnd + 5;
        }else{
          cwnd = cwnd + cwnd;
        }
      }

      float w = 1.0/16;
      double tmpbeta = 1.0 * col_count / awnd;
      beta = w * tmpbeta + (1-w) * beta;

      // double urgent = 1;
      float thre = 0.05;
      if (col_count > 0 && beta > thre){
        // awnd = awnd * (1 - beta/2);
        awnd = awnd * (1 - pow( (beta - thre)/(1-thre), urgent)/2);
      }else{
        awnd = awnd + 5;
      }


      awnd = awnd < awndmax_ ? awnd : awndmax_;
      awnd = awnd > awndmin_ ? awnd : awndmin_;

      
      awnd = awnd > cwnd ? cwnd : awnd; 

      col_count = 0;
      ecn_count = 0;
  }else{
    ;
  }

  cwnd = cwnd > 0 ? cwnd : 1;
  cwnd = cwnd < cwndmax_ ? cwnd : cwndmax_;
}
      


void simplecc::updatecwnd(int operate){
  if(enableasycc){
    cwnd += operate;
    cwnd = cwnd > 0 ? cwnd : 1;
    cwnd = cwnd < cwndmax_ ? cwnd : cwndmax_;
  }
}

uint32_t
simplecc::get_send_left_cound(){
      return last_send+1;
    }

uint32_t
simplecc::get_send_right_cound(){
      return lastack + 1 + cwnd;
    }

 void 
simplecc::update_ecn(uint8_t ecn){
      isecn = ecn;
    }

void 
simplecc::clean_ecn(){
      isecn = 0;
    }
void 
simplecc::setlast_send(uint32_t seq){
      last_send = seq;
}
void 
simplecc::setlast_ack(uint32_t ack){
      lastack = ack;
}
uint32_t 
simplecc::get_cwnd(){
  return cwnd;
}



} // Namespace ns3
