#!/usr/bin/env python3
import random
import os
import argparse

# 参数设置
parser = argparse.ArgumentParser(description='Generate topology files with different tail intensities')
parser.add_argument('--tail_intensity', type=str, default='low', choices=['low', 'medium', 'high'], help='Tail intensity level')
parser.add_argument('--output_dir', type=str, default='./lzy_mix', help='Output directory')
args = parser.parse_args()

# 确保输出目录存在
os.makedirs(args.output_dir, exist_ok=True)
os.makedirs(os.path.join(args.output_dir, "topology"), exist_ok=True)

# 为不同的拖尾强度设置参数
tail_params = {
    'low': {
        'bandwidth_variation': 0.1,  # 带宽变化范围为10%
        'delay_variation': 2,        # 延迟变化范围为基础延迟的2倍
        'affected_hosts_percent': 0.1, # 10%的主机有拖尾问题
    },
    'medium': {
        'bandwidth_variation': 0.3,  # 带宽变化范围为30%
        'delay_variation': 5,        # 延迟变化范围为基础延迟的5倍
        'affected_hosts_percent': 0.25, # 25%的主机有拖尾问题
    },
    'high': {
        'bandwidth_variation': 0.5,  # 带宽变化范围为50%
        'delay_variation': 10,       # 延迟变化范围为基础延迟的10倍
        'affected_hosts_percent': 0.4, # 40%的主机有拖尾问题
    }
}

# 获取当前拖尾强度的参数
tail_config = tail_params[args.tail_intensity]

# 大规模仿真的拓扑参数
NUM_TOR_SWITCHES = 8
NUM_CORE_SWITCHES = 8
NUM_HOSTS_PER_TOR = 24
TOTAL_NODES = NUM_TOR_SWITCHES + NUM_CORE_SWITCHES + NUM_TOR_SWITCHES * NUM_HOSTS_PER_TOR
FIRST_HOST_ID = NUM_TOR_SWITCHES + NUM_CORE_SWITCHES  # 第一个主机的ID是16

# 生成带拖尾的大规模拓扑文件
def generate_large_scale_topology():
    # 计算受拖尾影响的主机数量
    num_affected_hosts = int(NUM_TOR_SWITCHES * NUM_HOSTS_PER_TOR * tail_config['affected_hosts_percent'])
    
    # 随机选择受影响的主机ID（从第一个主机开始的索引）
    host_ids = list(range(FIRST_HOST_ID, TOTAL_NODES))
    affected_host_ids = random.sample(host_ids, num_affected_hosts)
    
    with open(f"{args.output_dir}/topology/topology_large_{args.tail_intensity}.txt", "w") as f:
        # 写入节点总数、交换机总数和链路总数
        f.write(f"# 拓扑：节点总数（{NUM_TOR_SWITCHES} 个 TOR 交换机 + {NUM_CORE_SWITCHES} 个核心交换机 + {NUM_TOR_SWITCHES} * {NUM_HOSTS_PER_TOR} 个节点 = {TOTAL_NODES} 个节点）\n")
        total_links = NUM_TOR_SWITCHES * NUM_CORE_SWITCHES + NUM_TOR_SWITCHES * NUM_HOSTS_PER_TOR
        f.write(f"{TOTAL_NODES} {NUM_TOR_SWITCHES + NUM_CORE_SWITCHES} {total_links}\n\n")
        
        # 写入节点编号说明
        f.write("# 节点编号说明:\n")
        f.write("# 0-7: TOR交换机\n")
        f.write("# 8-15: 核心交换机\n")
        for i in range(NUM_TOR_SWITCHES):
            start_id = FIRST_HOST_ID + i * NUM_HOSTS_PER_TOR
            end_id = start_id + NUM_HOSTS_PER_TOR - 1
            f.write(f"# {start_id}-{end_id}: 连接到TOR交换机{i}的主机\n")
        f.write("\n")
        
        # 写入交换机节点ID列表
        f.write("# 交换机节点ID列表\n")
        switch_ids = list(range(NUM_TOR_SWITCHES + NUM_CORE_SWITCHES))
        f.write(" ".join(map(str, switch_ids)) + "\n\n")
        
        # 写入TOR交换机连接到核心交换机的链路
        f.write("# TOR交换机连接到核心交换机\n")
        for tor in range(NUM_TOR_SWITCHES):
            for core in range(NUM_CORE_SWITCHES):
                core_id = NUM_TOR_SWITCHES + core
                f.write(f"{tor} {core_id} 100Gbps 8us 0\n")
            f.write("\n")
        
        # 写入主机连接到TOR交换机的链路
        for tor in range(NUM_TOR_SWITCHES):
            f.write(f"# 节点连接到 TOR 交换机 {tor}\n")
            for i in range(NUM_HOSTS_PER_TOR):
                host_id = FIRST_HOST_ID + tor * NUM_HOSTS_PER_TOR + i
                
                # 检查该主机是否受拖尾影响
                if host_id in affected_host_ids:
                    # 计算受影响的带宽和延迟
                    base_bw = 100
                    bw_variation = base_bw * tail_config['bandwidth_variation']
                    bandwidth = max(10, base_bw - random.uniform(0, bw_variation))  # 至少保留10Gbps
                    
                    base_delay = 4
                    delay_variation = base_delay * tail_config['delay_variation']
                    delay = base_delay + random.uniform(0, delay_variation)  # 增加延迟
                    
                    f.write(f"{tor} {host_id} {bandwidth:.0f}Gbps {delay:.0f}us 0\n")
                else:
                    # 正常带宽和延迟
                    f.write(f"{tor} {host_id} 100Gbps 4us 0\n")
            f.write("\n")

# 生成带拖尾的测试拓扑文件
def generate_test_topology():
    NUM_NODES = 40
    NUM_AFFECTED = int(NUM_NODES * tail_config['affected_hosts_percent'])
    
    # 随机选择受影响的主机ID
    host_ids = list(range(9, 9 + NUM_NODES))
    affected_host_ids = random.sample(host_ids, NUM_AFFECTED)
    
    with open(f"{args.output_dir}/topology/testtopo_{args.tail_intensity}.txt", "w") as f:
        # 写入节点总数、交换机总数和链路总数
        f.write(f"# 拓扑：测试版本（1 个 TOR 交换机 + 8 个核心交换机 + {NUM_NODES} 个节点）\n")
        total_nodes = 1 + 8 + NUM_NODES
        total_links = 8 + NUM_NODES  # 1个TOR连8个核心 + 40个主机连TOR
        f.write(f"{total_nodes} 9 {total_links}\n\n")
        
        # 写入节点编号说明
        f.write("# 节点编号说明:\n")
        f.write("# 0: TOR交换机\n")
        f.write("# 1-8: 核心交换机\n")
        f.write(f"# 9-{8+NUM_NODES}: 连接到TOR交换机的主机\n\n")
        
        # 写入交换机节点ID列表
        f.write("# 交换机节点ID列表\n")
        f.write("0 1 2 3 4 5 6 7 8\n\n")
        
        # 写入TOR交换机连接到核心交换机的链路
        f.write("# TOR交换机连接到核心交换机\n")
        for i in range(1, 9):
            f.write(f"0 {i} 100Gbps 2us 0\n")
        f.write("\n")
        
        # 写入受拖尾影响的节点
        f.write(f"# 节点连接到 TOR 交换机 0 - {NUM_AFFECTED}个有拖尾问题的节点\n")
        for host_id in sorted(affected_host_ids):
            # 计算受影响的带宽和延迟
            base_bw = 100
            bw_variation = base_bw * tail_config['bandwidth_variation']
            bandwidth = max(10, base_bw - random.uniform(0, bw_variation))  # 至少保留10Gbps
            
            base_delay = 2
            delay_variation = base_delay * tail_config['delay_variation']
            delay = base_delay + random.uniform(0, delay_variation)  # 增加延迟
            
            f.write(f"0 {host_id} {bandwidth:.0f}Gbps {delay:.0f}us 0\n")
        f.write("\n")
        
        # 写入正常节点
        f.write("# 其余节点正常带宽\n")
        for i in range(9, 9 + NUM_NODES):
            if i not in affected_host_ids:
                f.write(f"0 {i} 100Gbps 2us 0\n")

# 生成拓扑文件
print(f"正在生成拖尾强度为 {args.tail_intensity} 的拓扑文件...")
generate_large_scale_topology()
generate_test_topology()
print(f"拓扑文件已生成到 {args.output_dir}/topology 目录。") 