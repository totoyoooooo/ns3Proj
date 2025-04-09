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

 #ifndef UDP_PS_H
 #define UDP_PS_H
 
 #include "ns3/application.h"
 #include "ns3/event-id.h"
 #include "ns3/ptr.h"
 #include "ns3/address.h"
 #include "ns3/traced-callback.h"
 
 namespace ns3 {
 
 class Socket;
 class Packet;
 
 /**
  * \ingroup applications 
  * \defgroup udpecho UdpEcho
  */
 
 /**
  * \ingroup udpecho
  * \brief A Udp Echo server
  *
  * Every packet received is sent back.
  */
 class UdpPs : public Application 
 {
 public:
   /**
    * \brief Get the type ID.
    * \return the object TypeId
    */
   static TypeId GetTypeId (void);
   UdpPs ();
   virtual ~UdpPs ();
 
 protected:
   virtual void DoDispose (void);
 
 private:
 
   virtual void StartApplication (void);
   virtual void StopApplication (void);
 
   /**
    * \brief Handle a packet reception.
    *
    * This function is called by lower layers.
    *
    * \param socket the socket the packet was received to.
    */
   void HandleRead (Ptr<Socket> socket);
 
   uint16_t m_port; //!< Port on which we listen for incoming packets.
   Ptr<Socket> m_socket; //!< IPv4 Socket
   Ptr<Socket> m_socket6; //!< IPv6 Socket
   Address m_local; //!< local multicast address
   uint64_t m_totalReceived;   // 总接收包数（所有包的total累加）
   uint64_t m_totalMerged;     // 合并包数（isMerged=1的total累加）
   uint32_t mergetimes;
   uint32_t m_received;
   uint64_t mx;
   Time totalTime;
   std:: map<int ,int> mp[4];
 
   std:: map<int ,int> merge[4];
 
   /// Callbacks for tracing the packet Rx events
   TracedCallback<Ptr<const Packet> > m_rxTrace;
 
   /// Callbacks for tracing the packet Rx events, includes source and destination addresses
   TracedCallback<Ptr<const Packet>, const Address &, const Address &> m_rxTraceWithAddresses;
 
   bool aggregate_pkt(uint32_t key);
   std::map<uint32_t, uint32_t> key_to_bitmap_index;
   uint32_t total_worker;
   void  outputthrought();
   // for cc
   void update_ecn(uint32_t key,uint8_t ecn_mark);
   std::map<uint32_t, uint32_t> key_to_ecn_mark;
 uint32_t m_count_sent;
 };
 
 } // namespace ns3
 
 #endif /* UDP_PS_H */
 
 