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

#ifndef UDP_ECHO_WORKER_H
#define UDP_ECHO_WORKER_H

#include "ns3/application.h"
#include "ns3/event-id.h"
#include "ns3/ptr.h"
#include "ns3/ipv4-address.h"
#include "ns3/traced-callback.h"

namespace ns3 {

class Socket;
class Packet;
//for simple cc
//lyj add 
class simplecc{
  public:
    simplecc();
    ~simplecc();
    void initcwnd(uint32_t cwnd, uint32_t cwndmax, uint32_t onoffcc, uint32_t onoffasycc, uint32_t used_aggr);
    void adjustcwnd(uint32_t onofflzycc, uint32_t onoffcc);
    void adjustawnd(float tmpbeta); //useless in a2tp
    void updatecwnd(int operate); //lzy
    uint32_t get_send_left_cound();
    uint32_t get_send_right_cound();
    void update_ecn(uint8_t ecn);
    void clean_ecn();
    void setlast_send(uint32_t seq);
    void setlast_ack(uint32_t ack);
    uint32_t get_cwnd();
    uint32_t get_awnd(){
      return awnd;
    }
    float lastdelta; //lzy
    float weight; //lzy
    
    float alpha;
    uint32_t ecn_count;

    float beta;
    float urgent;
    uint32_t col_count;
    uint32_t recv_ack_num;
    
  private:
    uint32_t cwnd;
    uint32_t lastack;
    uint32_t last_send;
    uint32_t cwndmax_;
    uint32_t ssthre_;
    uint32_t awndssthre_;
    uint8_t isecn;
    uint32_t enablecc;
    uint32_t enableasycc; //lzy
    

    uint32_t awnd;
    uint32_t awndmax_;
    uint32_t awndmin_;




};


/**
 * \ingroup udpecho
 * \brief A Udp Echo client
 *
 * Every packet sent should be returned by the server and received here.
 */
class UdpWorker : public Application 
{
public:
  /**
   * \brief Get the type ID.
   * \return the object TypeId
   */
  static TypeId GetTypeId (void);

  UdpWorker ();

  virtual ~UdpWorker ();

  /**
   * \brief set the remote address and port
   * \param ip remote IP address
   * \param port remote port
   */
  void SetRemote (Address ip, uint16_t port);
  /**
   * \brief set the remote address
   * \param addr remote address
   */
  void SetRemote (Address addr);

  /**
   * Set the data size of the packet (the number of bytes that are sent as data
   * to the server).  The contents of the data are set to unspecified (don't
   * care) by this call.
   *
   * \warning If you have set the fill data for the echo client using one of the
   * SetFill calls, this will undo those effects.
   *
   * \param dataSize The size of the echo data you want to sent.
   */
  void SetDataSize (uint32_t dataSize);

  /**
   * Get the number of data bytes that will be sent to the server.
   *
   * \warning The number of bytes may be modified by calling any one of the 
   * SetFill methods.  If you have called SetFill, then the number of 
   * data bytes will correspond to the size of an initialized data buffer.
   * If you have not called a SetFill method, the number of data bytes will
   * correspond to the number of don't care bytes that will be sent.
   *
   * \returns The number of data bytes.
   */
  uint32_t GetDataSize (void) const;
  /**
   * Set the data fill of the packet (what is sent as data to the server) to 
   * the zero-terminated contents of the fill string string.
   *
   * \warning The size of resulting echo packets will be automatically adjusted
   * to reflect the size of the fill string -- this means that the PacketSize
   * attribute may be changed as a result of this call.
   *
   * \param fill The string to use as the actual echo data bytes.
   */
  void SetFill (std::string fill);

  /**
   * Set the data fill of the packet (what is sent as data to the server) to 
   * the repeated contents of the fill byte.  i.e., the fill byte will be 
   * used to initialize the contents of the data packet.
   * 
   * \warning The size of resulting echo packets will be automatically adjusted
   * to reflect the dataSize parameter -- this means that the PacketSize
   * attribute may be changed as a result of this call.
   *
   * \param fill The byte to be repeated in constructing the packet data..
   * \param dataSize The desired size of the resulting echo packet data.
   */
  void SetFill (uint8_t fill, uint32_t dataSize);

  /**
   * Set the data fill of the packet (what is sent as data to the server) to
   * the contents of the fill buffer, repeated as many times as is required.
   *
   * Initializing the packet to the contents of a provided single buffer is 
   * accomplished by setting the fillSize set to your desired dataSize
   * (and providing an appropriate buffer).
   *
   * \warning The size of resulting echo packets will be automatically adjusted
   * to reflect the dataSize parameter -- this means that the PacketSize
   * attribute of the Application may be changed as a result of this call.
   *
   * \param fill The fill pattern to use when constructing packets.
   * \param fillSize The number of bytes in the provided fill pattern.
   * \param dataSize The desired size of the final echo data.
   */
  void SetFill (uint8_t *fill, uint32_t fillSize, uint32_t dataSize);

  bool updateurgent;
  uint32_t m_fixedWindow; // 固定发送窗口大小

protected:
  virtual void DoDispose (void);

private:

  virtual void StartApplication (void);
  virtual void StopApplication (void);

  /**
   * \brief Schedule the next packet transmission
   * \param dt time interval between packets.
   */
  void ScheduleTransmit (Time dt);
  /**
   * \brief Send a packet
   */
  void Send (void);

  /**
   * \brief Handle a packet reception.
   *
   * This function is called by lower layers.
   *
   * \param socket the socket the packet was received to.
   */
  void HandleRead (Ptr<Socket> socket);

  uint32_t m_count; //!< Maximum number of packets the application will send
  Time m_interval; //!< Packet inter-send time
  uint32_t m_size; //!< Size of the sent packet

  uint32_t m_dataSize; //!< packet payload size (must be equal to m_size)
  uint8_t *m_data; //!< packet payload data

  uint32_t m_sent; //!< Counter for sent packets
  Ptr<Socket> m_socket; //!< Socket
  Address m_peerAddress; //!< Remote peer address
  uint16_t m_peerPort; //!< Remote peer port
  EventId m_sendEvent; //!< Event to send the next packet

  /// Callbacks for tracing the packet Tx events
  TracedCallback<Ptr<const Packet> > m_txTrace;

  /// Callbacks for tracing the packet Rx events
  TracedCallback<Ptr<const Packet> > m_rxTrace;
  
  /// Callbacks for tracing the packet Tx events, includes source and destination addresses
  TracedCallback<Ptr<const Packet>, const Address &, const Address &> m_txTraceWithAddresses;
  
  /// Callbacks for tracing the packet Rx events, includes source and destination addresses
  TracedCallback<Ptr<const Packet>, const Address &, const Address &> m_rxTraceWithAddresses;
  /// lyj add
  uint32_t m_cwnd; // send window
  uint32_t m_total_size; // Maximum bytes of the application will send
  uint32_t m_size_of_per_pkt; // Bytes of each pkt
  uint32_t m_max_index; //maximum index
  uint32_t m_complete; //count for aggreate pkt
  uint32_t m_count_sent;
  uint8_t m_host_id;
  uint16_t m_appid;
  uint8_t m_host_num;
  Time time;
  Time last_ack_time;
  uint32_t mx;

  void InitializeSending(void);
  void SwitchmlSend (void);
  void SendPacket (uint32_t index, uint8_t host_id, uint32_t bpAggr, uint32_t syn_urgent, uint8_t set_ecn);
  void SwitchmlHandleRead (Ptr<Socket> socket);
  void setindex (void);
  void logthroughput (void);
  void outputthrought(void);
  uint32_t UpdateNewLastACK(uint32_t ack_idx);

  std::map<uint32_t, bool> ack_bitmap;
  uint32_t last_cwnd_pkt;
  uint32_t last_awnd_pkt;
  uint32_t max_ack;
  uint32_t last_ack;
  uint32_t last_mack;

  uint32_t curseq_;
  simplecc m_cc;
  uint32_t onoffcc;
  uint32_t onoffasycc;
  uint32_t onofflzycc;
  uint32_t onoffawndcc;
  uint32_t used_aggr;
  uint32_t m_cwndmax;

  Time totalTime;
  double last_awnd_time;
  double last_urgent_time;
  double baseRTT; 

  }; // namespace ns3
}
#endif /* UDP_ECHO_WORKER_H */
