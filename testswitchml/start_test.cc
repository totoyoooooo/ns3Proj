/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
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

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/traffic-control-module.h"
#include "ns3/stats-module.h"
#include "ns3/flow-monitor-module.h"
#include <string>
#include <math.h>
#include <fstream>
#include <map>

// Default Network Topology
//
//       10.1.1.0
// n0 -------------- n1
//    point-to-point
//
// we orighin have 8*256=2048
// 1024 need to full 100Gbpss
//this is a test file
std::size_t numberworker = 4; // 2 jobs, and 2 workers of each job  
int numberaggregator = 2;
int pool_size = 100000;
int max_pool_size = 10000000;
int cwnd = 20;
int maxcwnd = 1000;
int ecn_thre = 256;
int onoffcc = 1;
int onoffasycc = 1;
int onofflzycc = 1;
int onoffawndcc = 1;
int onoffps = 0;
int onofftrace = 0;
int onofftimewindow = 0;
//int model_size = 528*1024*1024;
//int model_size = 528*1024*12;
int coloactednumber = 1;
uint32_t onsaatp = 0;
using namespace ns3;

std::string topofile = "./lzy_mix/topology_simple.txt";
std::string jobfile = "./lzy_mix/job_simple.txt";
std::string bgflowfile = "./lzy_mix/bgflow_simple.txt";
std::string trace_file = "./lzy_mix/lzytrace";

std::ifstream topof, jobf, bgf;

std::string max_buffer = "20000p";

// Function to sample flow statistics periodically
void SampleFlowStats(Ptr<FlowMonitor> monitor, std::ofstream* outFile, double interval)
{
  monitor->CheckForLostPackets();
  
  // Get current time
  double now = Simulator::Now().GetSeconds();
  
  // Calculate total traffic in this sampling interval
  std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();
  double totalThroughputMbps = 0.0;
  
  for (auto it = stats.begin(); it != stats.end(); ++it) {
    // Calculate instantaneous throughput based on received bytes
    if (it->second.rxBytes > 0) {
      // Calculate bytes received in last interval
      totalThroughputMbps += it->second.rxBytes * 8.0 / 1024 / 1024;
    }
  }
  
  // Write time and throughput to file
  *outFile << now << "\t" << totalThroughputMbps << std::endl;
  
  // Schedule next sampling
  Simulator::Schedule(Seconds(interval), &SampleFlowStats, monitor, outFile, interval);
}

// the simple cc can convergen to 8396Mbps when 
// cwnd = 100, poll_size = 1000, rates= "10Gbps", threshold = 4, Delay = 20us
// and the rat      e is 2592Mbps when using fixed cwnd
NS_LOG_COMPONENT_DEFINE ("FirstScriptExample");

int
main (int argc, char *argv[])
{
  CommandLine cmd (__FILE__);
  // cmd.AddValue ("low_link_speed", "straggler level", low_link_speed);
  // cmd.AddValue ("number_of_nodes", "number of nodes", numberworker);

  Time::SetResolution (Time::NS);
  // LogComponentEnable ("UdpWorkerApplication", LOG_LYJ);
  //LogComponentEnable ("UdpPsApplication", LOG_LYJ);
  // LogComponentEnable ("UdpWorkerApplication", LOG_LYJ);
  // LogComponentEnable ("UdpAggregatorApplication", LOG_LYJ);
  //LogComponentEnable("FifoQueueDisc",LOG_LYJ);
  //LogComponentEnable("Ipv4QueueDiscItem",LOG_LYJ);

  std::string configpath="./lzy_mix/config/testa2tp.txt";
  std::string jobpath="./lzy_mix/job/testjob.txt";
  int cmd_max_pool_size = 0;
  std::string modelType = "ResNet"; // Default model type
  std::string tailIntensity = "normal"; // Default tail intensity
  bool enableFlowMonitor = true; // Enable flow monitoring by default
  
  cmd.AddValue ("configpath", "config path", configpath);
  cmd.AddValue ("jobpath", "jobpath", jobpath);
  cmd.AddValue ("cmd_poolsize", "cmd_max_pool_size", cmd_max_pool_size);
  cmd.AddValue ("model", "Model type (ResNet or VGG)", modelType);
  cmd.AddValue ("tail", "Tail intensity (normal, medium, high)", tailIntensity);
  cmd.AddValue ("flowmon", "Enable flow monitoring", enableFlowMonitor);
  cmd.Parse (argc, argv);
  std::cout<<"config path "<<configpath<<std::endl;

  

  if (argc > 1) {
    std::ifstream conf;
    conf.open(configpath);
		while (!conf.eof()){
			std::string key;
			conf >> key;
			//std::cout << conf.cur << "\n";

			if (key.compare("ENABLE_CC") == 0){
				uint32_t v;
				conf >> v;
				onoffcc = v;
				if (onoffcc)
					std::cout << "ENABLE_CC\t" << "Yes" << "\n";
				else
					std::cout << "ENABLE_CC\t" << "No" << "\n";
			}else if (key.compare("ENABLE_ASYCC") == 0){
				uint32_t v;
				conf >> v;
				onoffasycc = v;
				if (onoffasycc)
					std::cout << "ENABLE_ASYCC\t" << "Yes" << "\n";
				else
					std::cout << "ENABLE_ASYCC\t" << "No" << "\n";
			}else if (key.compare("ENABLE_LZYCC") == 0){
				uint32_t v;
				conf >> v;
				onofflzycc = v;
				if (onofflzycc)
					std::cout << "ENABLE_LZYCC\t" << "Yes" << "\n";
				else
					std::cout << "ENABLE_LZYCC\t" << "No" << "\n";
			}else if (key.compare("ENABLE_AWNDCC") == 0){
				uint32_t v;
				conf >> v;
				onoffawndcc = v;
				if (onoffawndcc)
					std::cout << "ENABLE_AWNDCC\t" << "Yes" << "\n";
				else
					std::cout << "ENABLE_AWNDCC\t" << "No" << "\n";
			}else if (key.compare("ENABLE_PS") == 0){
				uint32_t v;
				conf >> v;
				onoffps = v;
				if (onoffps)
					std::cout << "ENABLE_PS\t" << "Yes" << "\n";
				else
					std::cout << "ENABLE_PS\t" << "No" << "\n";
			}else if (key.compare("INIT_CWND") == 0){
				uint32_t v;
				conf >> v;
				cwnd = v;
				std::cout << "INIT_CWND\t" << cwnd << "\n";
			}else if (key.compare("MAX_CWND") == 0){
				uint32_t v;
				conf >> v;
				maxcwnd = v;
				std::cout << "MAX_CWND\t" << maxcwnd << "\n";
			}else if (key.compare("MAX_BUFFER") == 0){
				std::string v;
				conf >> v;
				max_buffer = v;
				std::cout << "MAX_BUFFER\t" << max_buffer << "\n";
			}else if (key.compare("ECN_THRE") == 0){
				uint32_t v;
				conf >> v;
				ecn_thre = v;
				std::cout << "ECN_THRE\t" << ecn_thre << "\n";
			}else if (key.compare("MAX_AGGR_POOL") == 0){
				int v;
				conf >> v;
				max_pool_size = v;
				std::cout << "MAX_AGGR_POOL\t" << max_pool_size << "\n";
			}else if (key.compare("TOPO_FILE") == 0){
				std::string v;
				conf >> v;
				topofile = v;
				std::cout << "TOPO_FILE\t" << topofile << "\n";
			}else if (key.compare("JOB_FILE") == 0){
				std::string v;
				conf >> v;
				jobfile = v;
				std::cout << "JOB_FILE\t" << jobfile << "\n";
			}else if (key.compare("BGFLOW_FILE") == 0){
				std::string v;
				conf >> v;
				bgflowfile = v;
				std::cout << "BGFLOW_FILE\t" << bgflowfile << "\n";
			}else if (key.compare("TRACE_FILE") == 0){
				std::string v;
				conf >> v;
				trace_file = v;
				std::cout << "TRACE_FILE\t" << trace_file << "\n";
      }else if (key.compare("ENABLE_TRACE") == 0){
				int v;
				conf >> v;
				onofftrace = v;
				std::cout << "ENABLE_TRACE\t" << onofftrace << "\n";
      } else if (key.compare("ENABLE_TIMEWINDOW") == 0){
        int v;
				conf >> v;
				onofftimewindow = v;
				std::cout << "ENABLE_TIMEWINDOW\t" << onofftimewindow << "\n";
      }
      
    }
  }
 
  
  //generate topo
  topof.open(topofile.c_str());

  uint32_t node_num, switch_num, link_num;
	topof >> node_num >> switch_num >> link_num;
  
  NodeContainer nodes;
  nodes.Create (node_num);

  std::vector<uint32_t> node_type(node_num, 0);

  struct Node_To_Node{
    uint32_t src;
    uint32_t dst;
    NetDeviceContainer dev;
    Ipv4InterfaceContainer ip;
  };
  std::vector<Node_To_Node> node_to_node;
  node_to_node.reserve(link_num);

  std::vector<std::vector<int> > node_link_table(node_num, std::vector<int>(node_num, -1));
  std::vector<ns3::Ipv4Address> node_ip(node_num);

	for (uint32_t i = 0; i < switch_num; i++){
		uint32_t sid;
		topof >> sid;
		node_type[sid] = 1;
	}

 // std::vector<NetDeviceContainer> node_to_node_dev;


  PointToPointHelper pointToPoint;

  for (uint32_t i = 0; i < link_num; i++){

    Node_To_Node temp;
    node_to_node.push_back(temp);


    uint32_t src, dst;
    std::string data_rate, link_delay;
    double error_rate;
    topof >> src >> dst >> data_rate >> link_delay >> error_rate;
          
    node_link_table[src][dst] = i;
    node_link_table[dst][src] = i;

    node_to_node[i].src = src;
    node_to_node[i].dst = dst;

    Ptr<Node> snode = nodes.Get(src), dnode = nodes.Get(dst);
    pointToPoint.SetDeviceAttribute ("DataRate", StringValue (data_rate));
    pointToPoint.SetChannelAttribute ("Delay", StringValue (link_delay));
    // pointToPoint.SetChannelAttribute ("Delay", StringValue ("2us"));

    node_to_node[i].dev = pointToPoint.Install (snode, dnode);
    // std::cout<<"snodeID "<<snode->GetId()<<" sDeviceID "\
    //   <<node_to_node[i].dev.Get(0)->GetIfIndex()\
    //   <<" dnodeID "<<dnode->GetId()<<" dDeviceID "\
    //   <<node_to_node[i].dev.Get(1)->GetIfIndex()<<std::endl;

    fflush(stdout);
  }

  InternetStackHelper stack;
  stack.InstallAll();

  Ipv4AddressHelper address;
  address.SetBase ("10.1.1.0", "255.255.255.0");
  address.NewNetwork ();
 // std::vector<Ipv4InterfaceContainer> node_to_node_dev_ip;
 // node_to_node_dev_ip.reserve(node_num);


  TrafficControlHelper tch;
  tch.SetRootQueueDisc ("ns3::FifoQueueDisc","MaxSize", StringValue (max_buffer),
                        "ECNThreshold", UintegerValue (ecn_thre));

  //generate ecn queue
  for (int i = 0; i < link_num; i++){

    int src = node_to_node[i].src;
    int dst = node_to_node[i].dst;

    // if (src == 264 && node_type[dst] == 0 && dst < 64 || src == 264 && node_type[dst] == 1 || src == 264 && node_type[dst] == 0 && dst < 64 ){
    //   ;
    // }

    if(node_type[src] == 1) tch.Install (node_to_node[i].dev.Get(0));
    if(node_type[dst] == 1) tch.Install (node_to_node[i].dev.Get(1));

    // if(node_type[src] == 1 && (dst>=32 && dst<=35 || dst == 39) && src!=39) tch.Install (node_to_node[i].dev.Get(0));
    // if(node_type[dst] == 1 && (src>=32 && src<=35 && dst!= 39 || src == 39)) tch.Install (node_to_node[i].dev.Get(1));


    node_to_node[i].ip = address.Assign(node_to_node[i].dev);
    if(node_type[src] == 0) node_ip[src] = node_to_node[i].ip.GetAddress(0);
    if(node_type[dst] == 0) node_ip[dst] = node_to_node[i].ip.GetAddress(1);
    //std::cout<<"node_to_node ["<<node_to_node[i].src<<","<<node_to_node[i].dst<< "] "<<node_to_node[i].ip.GetAddress(0)<<" <--> "<< node_to_node[i].ip.GetAddress(1)<<std::endl;
    //std::cout<<"mac node_to_node ["<<node_to_node[i].src<<","<<node_to_node[i].dst<< "] "<<node_to_node[i].dev.Get(0)<<" <--> "<< node_to_node[i].dev.Get(1)<<std::endl;
    address.NewNetwork ();

  }
  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();
  std::string probeType;
  std::string tracePath;
  FileHelper fileHelper;
  if (onofftrace){

    probeType = "ns3::PacketProbe";
    tracePath = "NodeList/*/DeviceList/*/$ns3::PointToPointNetDevice/TxQueue/Dequeue";


    // Configure the file to be written, and the formatting of output data.
    fileHelper.ConfigureFile ((const std::string)trace_file,
                              FileAggregator::FORMATTED);

    // Set the labels for this formatted output file.
    fileHelper.Set2dFormat ("Time(Seconds)= %.4e\tPacketByteCount= %.0f");

    // Specify the probe type, trace source path (in configuration namespace), and
    // probe output trace source ("OutputBytes") to write.
    fileHelper.WriteProbe (probeType,
                          tracePath,
                          "OutputBytes");
  }
  
  //ipv4_global_routing.cc enable ecmp
  // AsciiTraceHelperForDevice ascihelper;
  // ascihelper.EnableAsciiIpv4 ("prefix", 21, 1);

  
  //generate job
  // jobf.open(jobfile.c_str());jobpath
  jobf.open(jobpath.c_str());
  uint32_t app_num;
	jobf >> app_num;
  std::vector<uint32_t> is_node_aggr(node_num, 0);
  uint16_t udpport = 9;
  ApplicationContainer clientApps[1024][64];
  ApplicationContainer serverApps[1024];
  Ptr<UdpAggregator> udpAggregator[1024];
  uint16_t psudport[1024];
  memset(psudport, 0 , sizeof(psudport));
  if (cmd_max_pool_size > 0){
    max_pool_size = cmd_max_pool_size;
  }
  for (int i=0; i < app_num; i++){

    // udpport++;
    // std::cout<<"udpport "<< udpport<<std::endl;
    uint32_t appid, workernum, aggregator_node, ps_node, max_count, aggr_used, model_size, cross_rack;
    float interval, start_time;
    jobf >> appid >> workernum >> aggregator_node >> ps_node >> max_count >> interval >> aggr_used >> model_size >> start_time >> cross_rack;
    printf("StartAPP %d start %lf \n", appid, start_time);
    // std::cout<<appid <<" "<< workernum <<" "<< aggregator_node <<" "<< ps_node <<" "<< max_count <<" "<< interval <<" "<< aggr_used <<" "<< model_size <<" "<< start_time<<std::endl;
    if (psudport[ps_node] == 0){
      psudport[ps_node] = udpport;
    }else{
      psudport[ps_node] ++;
    }
    if(is_node_aggr[aggregator_node] == 0){
      is_node_aggr[aggregator_node] = 1;
      UdpAggregatorHelper aggregator (udpport);
      aggregator.SetAttribute ("Msize",UintegerValue(max_pool_size));
      aggregator.SetAttribute ("MAX_Host_Number",UintegerValue(64));
      //aggregator.SetAttribute ("Aggregateid",UintegerValue(0));
      //aggregator.SetAttribute ("RemoteAddress",AddressValue(aggregator_to_sink_ip.GetAddress(1)));
      aggregator.SetAttribute ("RemotePort",UintegerValue(udpport));
      aggregator.SetAttribute ("Port",UintegerValue(udpport));
      aggregator.SetAttribute ("Level",UintegerValue(1));
      aggregator.SetAttribute("OnSaatp",UintegerValue(onsaatp));
      aggregator.SetAttribute("ONOFFASYCC",UintegerValue(onoffasycc));
      aggregator.SetAttribute("ONOFFPS",UintegerValue(onoffps));
      aggregator.SetAttribute("ONOFFTIMEWINDOW",UintegerValue(onofftimewindow));
      aggregator.SetAttribute("TimeWindow",DoubleValue(0.000001));
       
      serverApps[aggregator_node] = aggregator.Install (nodes.Get(aggregator_node));
      udpAggregator[aggregator_node] = DynamicCast<UdpAggregator> (serverApps[aggregator_node].Get(0));
      serverApps[aggregator_node].Start (Seconds (1.0));
      serverApps[aggregator_node].Stop (Seconds (100.0)); 
    }
    std::vector<ns3::Address> remotes;
    std::vector<uint16_t> port; 
    // for(int j=0; j < node_num; j++){
    //   int link_index = node_link_table[aggregator_node][j];
    //   if(link_index >= 0 && node_type[j] == 0){
    //     //int hostip_index = node_to_node[link_index].src == aggregator_node ? 1 : 0;
    //     //remotes.push_back(InetSocketAddress(node_to_node[link_index].ip.GetAddress(hostip_index),9));
    //     remotes.push_back(InetSocketAddress(node_ip[j],udpport));
    //     port.push_back(udpport);
    //     //std::cout<<node_to_node[link_index].ip.GetAddress(1)).GetIpv4()<<std::endl;
    //   }

    // }
    int link_index = node_link_table[aggregator_node][ps_node];
    if(link_index >= 0 && node_type[ps_node] == 0){
      // printf("");
      remotes.push_back(InetSocketAddress(node_ip[ps_node],psudport[ps_node]));
      port.push_back(psudport[ps_node]);
    }else{
      printf("error ps_node\n");
    }
    printf("debug setremotes remotesize %d port size %d appid %d\n", remotes.size(), port.size(), appid);
    udpAggregator[aggregator_node]->setremotes(remotes,port,appid);  

    
    UdpPsHelper PSer (psudport[ps_node]);
    PSer.SetAttribute ("Toalworker",UintegerValue(workernum));
    ApplicationContainer PsApps = PSer.Install (nodes.Get(ps_node));
    PsApps.Start (Seconds (1.0));
    PsApps.Stop (Seconds (100.0));

    ns3::Ipv4Address aggr_address;
    link_index = node_link_table[aggregator_node][ps_node];
    int portip_index = node_to_node[link_index].src == aggregator_node ? 0 : 1;
    aggr_address = node_to_node[link_index].ip.GetAddress(portip_index);
    for(std::size_t j = 0; j < workernum; j++){
      int accwnd = maxcwnd;
      int acpool_size = max_pool_size;
      double basertt = 0.000016;
      if (!cross_rack){
        accwnd = accwnd/2;
        basertt = basertt/2;
      }

      if (!onoffps){
        int allc = max_pool_size  / (app_num);
        if(!cross_rack){
          accwnd = maxcwnd / 4;
          basertt = 0.000004;
        }else{
          accwnd = maxcwnd / 4 * 3;
          basertt = 0.000012;
        }
        accwnd = accwnd > allc / 2 ? allc / 2 : accwnd;
        acpool_size = accwnd;
      }
      int worker_node;
      jobf >> worker_node;
      UdpWorkerHelper workers (aggr_address, udpport);
      workers.SetAttribute ("MaxPackets", UintegerValue (max_count));
      workers.SetAttribute ("Interval", TimeValue (Seconds (interval)));
      workers.SetAttribute ("PacketSize", UintegerValue (256));
      // workers.SetAttribute ("Maxbytes", UintegerValue (model_size*1000000));
      workers.SetAttribute ("Maxbytes", UintegerValue (model_size*10000));
      workers.SetAttribute ("CWND", UintegerValue (accwnd));
      workers.SetAttribute ("MAXCWND", UintegerValue (accwnd));
      workers.SetAttribute ("Host_number",UintegerValue(workernum));
      workers.SetAttribute ("HOSTID", UintegerValue (j));
      workers.SetAttribute ("APPID", UintegerValue (appid));
      workers.SetAttribute ("ONOFFCC", UintegerValue (onoffcc));
      workers.SetAttribute ("ONOFFAWNDCC", UintegerValue(onoffawndcc));
      workers.SetAttribute ("USEDAGGR", UintegerValue (acpool_size));
      workers.SetAttribute ("BASERTT", DoubleValue(basertt));
      
      clientApps[appid][j] = workers.Install (nodes.Get (worker_node));
      std::cout<<"app"<<appid<<" wnode "<<worker_node <<"-->"<<aggregator_node<<" ip "<<aggr_address<<std::endl;
      clientApps[appid][j].Start (Seconds (start_time));
      clientApps[appid][j].Stop (Seconds (100.0)); 
    }


    
    fflush(stdout);
  }

  //generate bgflow
//   #bg_num
// #[nouse0, nouse1, src_node, dst_node, pkt_count, start_time]


  uint32_t dlport = 9999;
  float k = 0;
  float interval = 0.024;
  float max_time = 0;
  while(k*interval <= max_time){
    bgf.open(bgflowfile.c_str());
    
    uint32_t bg_num;
    bgf >> bg_num;
    //std::cout<<bg_num<<std::endl;
    for( int i = 0; i < bg_num; i++){
  
      uint32_t nouse0, nouse1, bgsrc, bgdst, pkt_num=0;
      dlport++;
      float start_time=0;
      bgf >> nouse0 >> nouse1 >> bgsrc >> bgdst >> pkt_num >> start_time;
      //std::cout<<nouse0<<" "<<nouse1<<" "<<bgsrc<<" "<<bgdst<<" "<<pkt_num<<" "<<start_time<<std::endl;
      // std::cout<<"time "<<start_time<<std::endl;
      start_time = start_time + k*interval;
      
      if(node_type[bgsrc] == 1 || node_type[bgdst] == 1) { std::cout<<"error the node is not a host in bgflow "<<i<<std::endl; }




      // Create an optional packet sink to receive these packets

      PacketSinkHelper server ("ns3::UdpSocketFactory", Address (InetSocketAddress (Ipv4Address::GetAny (), dlport)));
      // PacketSinkHelper server ("ns3::TcpSocketFactory", Address (InetSocketAddress (Ipv4Address::GetAny (), dlport)));
      
      // UdpServerHelper server (dlport);
      ApplicationContainer udpserverapps = server.Install (nodes.Get(bgdst));
      udpserverapps.Start (Seconds (1.0));
      udpserverapps.Stop (Seconds (100.0));

      // std::cout<<"node "<<bgdst<<" ip "<<node_ip[bgdst]<<std::endl;

      // BulkSendHelper client ("ns3::TcpSocketFactory", InetSocketAddress (node_ip[bgdst], dlport));
      // client.SetAttribute ("MaxBytes", UintegerValue (pkt_num));
      OnOffHelper onoff ("ns3::UdpSocketFactory", 
                        Address (InetSocketAddress (node_ip[bgdst], dlport)));
      onoff.SetConstantRate (DataRate ("100Gib/s"), 1024);

      ApplicationContainer udpclientapps = onoff.Install (nodes.Get(bgsrc));
      // Start the application
      pkt_num = pkt_num/1024 + 1;
      // onoff.SetAttribute ("MaxPackets", UintegerValue (pkt_num));
      // onoff.SetAttribute ("Interval", TimeValue (Seconds(0)));
      onoff.SetAttribute ("PacketSize", UintegerValue (1024));
      udpclientapps.Start (Seconds (start_time));
      udpclientapps.Stop (Seconds (100.0));


      fflush(stdout);
    }
    bgf.close();
    k++;    
  }

  Simulator::Stop(Seconds(100.1));

  // After Simulator::Stop but before Simulator::Run, set up flow monitoring
  // Initialize Flow Monitor
  FlowMonitorHelper flowHelper;
  Ptr<FlowMonitor> flowMonitor;
  std::ofstream flowFile;
  std::map<std::string, std::ofstream> timeSeriesFiles;
  
  if (enableFlowMonitor) {
    flowMonitor = flowHelper.InstallAll();
    
    // Create the main flow stats file
    flowFile.open("Flow.txt", std::ios::out);
    flowFile << "# Flow statistics for model: " << modelType << ", tail intensity: " << tailIntensity << std::endl;
    flowFile << "# FlowID\tSrcIP\tDstIP\tTxPackets\tRxPackets\tLostPackets\tDelaySum(s)\tJitterSum(s)\tThroughput(Mbps)" << std::endl;
    
    // Create time series file for this model and tail intensity
    std::string timeSeriesFileName = modelType + "_" + tailIntensity + "_traffic.txt";
    timeSeriesFiles[timeSeriesFileName].open(timeSeriesFileName, std::ios::out);
    timeSeriesFiles[timeSeriesFileName] << "# Time(s)\tTraffic(Mbps)" << std::endl;
    
    // Set up periodic flow sampling for time series data
    Simulator::Schedule(Seconds(0.01), &SampleFlowStats, flowMonitor, &timeSeriesFiles[timeSeriesFileName], 0.01);
  }

  Simulator::Run ();
  
  // After simulation finishes, collect and save flow statistics
  if (enableFlowMonitor) {
    flowMonitor->CheckForLostPackets();
    
    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier>(flowHelper.GetClassifier());
    std::map<FlowId, FlowMonitor::FlowStats> stats = flowMonitor->GetFlowStats();
    
    for (auto it = stats.begin(); it != stats.end(); ++it) {
      Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow(it->first);
      
      // Calculate metrics
      FlowId flowId = it->first;
      uint32_t txPackets = it->second.txPackets;
      uint32_t rxPackets = it->second.rxPackets;
      uint32_t lostPackets = txPackets - rxPackets;
      double delaySum = it->second.delaySum.GetSeconds();
      double jitterSum = it->second.jitterSum.GetSeconds();
      
      double throughput = 0;
      if (it->second.rxBytes > 0 && 
          it->second.timeLastRxPacket.GetSeconds() > it->second.timeFirstTxPacket.GetSeconds()) {
        throughput = it->second.rxBytes * 8.0 / 
                    (it->second.timeLastRxPacket.GetSeconds() - it->second.timeFirstTxPacket.GetSeconds()) / 1024 / 1024;
      }
      
      // Write to file
      flowFile << flowId << "\t"
               << t.sourceAddress << "\t"
               << t.destinationAddress << "\t"
               << txPackets << "\t"
               << rxPackets << "\t"
               << lostPackets << "\t"
               << delaySum << "\t"
               << jitterSum << "\t"
               << throughput << std::endl;
    }
    
    flowFile.close();
    
    // Close time series files
    for (auto& file : timeSeriesFiles) {
      file.second.close();
    }
    
    // Generate a summary file with model and tail intensity information
    std::ofstream summaryFile("FlowSummary.txt", std::ios::app);
    summaryFile << "Model: " << modelType << ", Tail Intensity: " << tailIntensity 
               << ", Total Flows: " << stats.size() << std::endl;
    summaryFile.close();
  }

  Simulator::Destroy ();
  return 0;
}

//ECMP not used
// std::vector<ns3::Ipv4Address> ecmpaddr;
// ecmpaddr.reserve(node_num);
// for(int j=0; j < node_num; j++){
//     int link_index = node_link_table[aggregator_node][j];
          
//     if(link_index >= 0 && node_type[j] == 1){
//       //std::cout << "find a path "<<node_to_node[link_index].src<<" "<<node_to_node[link_index].dst<<std::endl;
//       int portip_index = node_to_node[link_index].src == aggregator_node ? 0 : 1;
//       ecmpaddr.push_back(node_to_node[link_index].ip.GetAddress(portip_index));
//       std::cout<<"ECMP "<<j<<"-->"<< aggregator_node <<" ip "<<node_to_node[link_index].ip.GetAddress(portip_index)<<std::endl;
//     }
// }

// int ecmpsize = ecmpaddr.size();
// if(workernum > ecmpsize){
//   std::cout<<"error workernum > ecmpsize "<<workernum<<" "<<ecmpsize<<std::endl;
// }