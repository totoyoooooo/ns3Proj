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

 #ifndef UDP_AGGREGATOR_H
 #define UDP_AGGREGATOR_H
 
 #include "ns3/application.h"
 #include "ns3/event-id.h"
 #include "ns3/ptr.h"
 #include "ns3/address.h"
 #include "ns3/traced-callback.h"
 
 #include <map>
 #include <vector>
 #include <fstream>
 #include <set>
 #include <queue>
 
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
 class UdpAggregator : public Application 
 {
 public:
   /**
    * \brief Get the type ID.
    * \return the object TypeId
    */
   static TypeId GetTypeId (void);
   UdpAggregator ();
   virtual ~UdpAggregator ();
 
   // void setremotes(std::vector<Address> address, std::vector<uint16_t> port);
   void setremotes(std::vector<Address> address, std::vector<uint16_t> port, uint32_t appid); 
 
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
 
   /// Callbacks for tracing the packet Rx events
   TracedCallback<Ptr<const Packet> > m_rxTrace;
 
   /// Callbacks for tracing the packet Rx events, includes source and destination addresses
   TracedCallback<Ptr<const Packet>, const Address &, const Address &> m_rxTraceWithAddresses;
 
 
   //lyj add
   uint32_t pool_size; // number of aggregator
   uint16_t max_host_num; //number of host
   uint32_t switch_level; 
   uint32_t aggregate_id;
   uint32_t minspace;
   Address m_peerAddress; //!< Remote peer address
   uint16_t m_peerPort; //!< Remote peer port
   // add
   struct AggregatorUnit {
    EventId timer;        // 定时器事件
    Time startTime;       // 聚合开始时间
    bool isFlushing = false; // 防重入标志
    Ptr<Packet> mergedPacket; // 合并后的数据包
  };
  std::map<uint16_t, std::map<uint32_t, AggregatorUnit>> m_units; 
  double m_timeWindow;    // 聚合时间窗口
  EventId m_sendEvent;         // 定时发送事件
  uint32_t m_merged_packet; //合并包数
  uint32_t m_merged; // 当前交换机合并次数
  uint32_t m_forwarded;  // 当前交换机转发次数
  uint32_t m_timeout_forwarded;  // 当前交换机转发次数

  std::map<std::pair<uint16_t, uint32_t>, std::vector<Time>> m_arrivalTimes;
  std::map<uint16_t, double> m_currentT;  // 各app的当前T值
  std::map<uint16_t, double> alpha;
  std::map<uint16_t, double> mu;
  uint32_t m_sampleThreshold;      // 触发参数估计的样本阈值
  double lambda;          // 基础到达率
  int Y = 250;    // 聚合器数量
  int n = 10;    // 组大小
  std::ofstream m_timeLogFile;      // 时间记录文件
  std::string m_timeDataFile; // 数据文件路径
  std::map<std::pair<uint16_t, uint32_t>, std::vector<double>> m_cachedSamples; // 缓存样本
  void EstimateAndUpdateT(uint16_t appid, uint32_t key); // 时间记录
  
  void LoadCachedSamples();// 样本加载
  std::vector<double> estimateLomaxParameters(const std::vector<double>& samples);
  void LogArrivalTime(uint16_t appId, uint32_t key);
  std::map<uint16_t, std::deque<Time>> m_arrivalWindows;  // 每个应用的到达时间窗口
  std::map<uint16_t, double> m_lambda;  // 每个应用的估计到达率
  uint32_t m_windowSize;  // 滑动窗口大小
  double m_ewmaAlpha;  // EWMA平滑因子
  std::map<uint16_t, std::set<uint8_t>> m_observedHosts;  // 每个应用观察到的主机ID
  void UpdateLambdaEstimate(uint16_t appId, Time arrivalTime);
   
   // std::vector< Ptr<Socket> > m_socket_for_up; //!< IPv4 Socket
   // std::vector< Address > multiple_address_for_up;
    // std::vector< uint16_t > multiple_port_for_up;

   std::map<uint32_t, Ptr<Socket> > m_socket_for_up;
   std::map<uint32_t, Address> multiple_address_for_up;
   std::map<uint32_t, uint16_t > multiple_port_for_up;
 
   std::map<uint16_t, std::map<uint32_t,uint32_t> > app_and_key_to_bitmap_index;
 
   std::map<uint16_t, std::map<uint32_t,uint32_t> > app_and_key_to_kicked; 
   std::map<uint16_t, std::map<uint32_t,uint32_t> > app_and_key_forwarded; 
 
   std::vector<uint32_t> unused;
 
   int64_t *start_time_index;
   uint32_t *delta_sum;
   uint32_t **delta_host; //lzy
   std::set<Address> group[100];
   bool **bitmap; 
   bool *aggr_collision; //lzy
   uint16_t *idx_appid; //lzy
   uint32_t *count_pkt;
   int64_t **interval;
 
   void initmem (uint32_t size, uint16_t host);
   void ForceFlush(uint16_t appid, uint32_t key);
   bool aggregate_pkt (uint16_t appid,uint32_t key, uint8_t hostid, uint8_t hostnum, Ptr<Packet> packet, Ptr<Socket> socket, Address from, uint8_t is_bpAggr, uint32_t syn_urgent);
   //bool aggregate_pkt (uint16_t appid,uint32_t key, uint32_t host,Ptr<Packet> packet);
   void broadcast(Ptr<Socket> socket,Ptr<Packet> packet,uint16_t appid,uint32_t key, uint8_t is_bpAggr, uint16_t is_Col, uint32_t syn_urgent, uint8_t set_ecn);
   void upload(Ptr<Packet> packet,uint16_t appid,uint32_t key, uint8_t is_bpAggr, uint32_t syn_urgent);
 
   void HandleReadForUp (Ptr<Socket> socket);
   void cleanaggregator(uint16_t appid,uint32_t key);
 
   std::vector<uint32_t> m_count_sent;
   std::vector<uint32_t> m_count_key;
   void outputthrought();
   uint32_t minimumkey;
   uint32_t maximumkey;
   uint32_t endkey;
   uint32_t on_saatp;
   uint32_t onoffasycc;
   uint32_t onoffps;
   uint32_t onofftimewindow;
   bool first;
   Time start_time;
 
   uint32_t collision_marking;
 
   bool hasForwarded(uint16_t appid, uint32_t key);
   bool iskicked(uint16_t appid, uint32_t key);
   bool isexist(uint16_t appid, uint32_t key);
   bool updateindexmap(uint16_t appid, uint32_t key, uint32_t index);
   
 };
 
 } // namespace ns3
 
 #endif /* UDP_AGGREGATOR_H */
 
 