#!/usr/bin/env python3
import random
import os
import argparse

# 参数设置
parser = argparse.ArgumentParser(description='Generate job configuration for large-scale simulation')
parser.add_argument('--jobs', type=int, default=20, help='Number of jobs to generate')
parser.add_argument('--tail_intensity', type=str, default='low', choices=['low', 'medium', 'high'], help='Tail intensity level')
parser.add_argument('--output_dir', type=str, default='./lzy_mix/job', help='Output directory')
parser.add_argument('--config_dir', type=str, default='./lzy_mix/config', help='Config directory')
args = parser.parse_args()

# 确保输出目录存在
os.makedirs(args.output_dir, exist_ok=True)
os.makedirs(args.config_dir, exist_ok=True)

# 配置参数
NUM_JOBS = args.jobs
HOST_START = 16  # 主机节点从16开始
HOST_END = 207   # 主机节点到207结束
NUM_HOSTS = HOST_END - HOST_START + 1  # 总共192个主机
MAX_JOBS_PER_PS = 8  # 每个PS最多支持8个job
TOR_SWITCH_START = 0  # TOR交换机从0开始
TOR_SWITCH_END = 7    # TOR交换机到7结束
HOSTS_PER_TOR = 24    # 每个TOR连接24个主机

# 为不同的拖尾强度设置参数
tail_params = {
    'low': {
        'bandwidth_variation': 0.1,  # 带宽变化范围为10%
        'interval_variation': 0.1,    # 发送间隔变化范围为10%
    },
    'medium': {
        'bandwidth_variation': 0.3,  # 带宽变化范围为30%
        'interval_variation': 0.3,    # 发送间隔变化范围为30%
    },
    'high': {
        'bandwidth_variation': 0.5,  # 带宽变化范围为50%
        'interval_variation': 0.5,    # 发送间隔变化范围为50%
    }
}

# 获取当前拖尾强度的参数
tail_config = tail_params[args.tail_intensity]

# 跟踪已分配的资源
ps_job_count = [0] * NUM_HOSTS  # 每个PS主机上运行的job数量
worker_allocation = [False] * NUM_HOSTS  # 每个主机是否已被分配为worker

# 创建配置文件 - ATP
with open(f"{args.config_dir}/testatp.txt", "w") as f:
    f.write("ENABLE_CC 1\n")
    f.write("ENABLE_ASYCC 0\n")
    f.write("ENABLE_LZYCC 0\n")
    f.write("ENABLE_AWNDCC 0\n")
    f.write("ENABLE_PS 1\n")
    f.write("INIT_CWND 340\n")
    f.write("MAX_CWND 340\n")
    f.write("MAX_BUFFER 50000p\n")
    f.write("ECN_THRE 4\n")
    f.write("MAX_AGGR_POOL 450\n")
    f.write(f"TOPO_FILE ./lzy_mix/topology/topology_large_{args.tail_intensity}.txt\n")
    f.write(f"JOB_FILE {args.output_dir}/large_jobs_{args.tail_intensity}.txt\n")
    f.write("BGFLOW_FILE ./lzy_mix/bgflow_large.txt\n")
    f.write("TRACE_FILE ./lzy_mix/lzytrace\n")
    f.write("ENABLE_TRACE 0\n")
    f.write("ENABLE_TIMEWINDOW 0\n")

# 创建配置文件 - A2TP
with open(f"{args.config_dir}/testa2tp.txt", "w") as f:
    f.write("ENABLE_CC 1\n")
    f.write("ENABLE_ASYCC 1\n")
    f.write("ENABLE_LZYCC 0\n")
    f.write("ENABLE_AWNDCC 0\n")
    f.write("ENABLE_PS 1\n")
    f.write("INIT_CWND 340\n")
    f.write("MAX_CWND 340\n")
    f.write("MAX_BUFFER 50000p\n")
    f.write("ECN_THRE 4\n")
    f.write("MAX_AGGR_POOL 450\n")
    f.write(f"TOPO_FILE ./lzy_mix/topology/topology_large_{args.tail_intensity}.txt\n")
    f.write(f"JOB_FILE {args.output_dir}/large_jobs_{args.tail_intensity}.txt\n")
    f.write("BGFLOW_FILE ./lzy_mix/bgflow_large.txt\n")
    f.write("TRACE_FILE ./lzy_mix/lzytrace\n")
    f.write("ENABLE_TRACE 0\n")
    f.write("ENABLE_TIMEWINDOW 0\n")

# 创建配置文件 - TimeWindow
with open(f"{args.config_dir}/testtime.txt", "w") as f:
    f.write("ENABLE_CC 1\n")
    f.write("ENABLE_ASYCC 0\n")
    f.write("ENABLE_LZYCC 0\n")
    f.write("ENABLE_AWNDCC 0\n")
    f.write("ENABLE_PS 1\n")
    f.write("INIT_CWND 340\n")
    f.write("MAX_CWND 340\n")
    f.write("MAX_BUFFER 50000p\n")
    f.write("ECN_THRE 4\n")
    f.write("MAX_AGGR_POOL 450\n")
    f.write(f"TOPO_FILE ./lzy_mix/topology/topology_large_{args.tail_intensity}.txt\n")
    f.write(f"JOB_FILE {args.output_dir}/large_jobs_{args.tail_intensity}.txt\n")
    f.write("BGFLOW_FILE ./lzy_mix/bgflow_large.txt\n")
    f.write("TRACE_FILE ./lzy_mix/lzytrace\n")
    f.write("ENABLE_TRACE 0\n")
    f.write("ENABLE_TIMEWINDOW 1\n")

# 辅助函数：获取主机在整个网络中的实际索引
def get_host_real_index(host_idx):
    """将0-191的主机索引转换为16-207的实际节点ID"""
    return HOST_START + host_idx

# 辅助函数：获取主机连接的TOR交换机ID
def get_host_tor(host_idx):
    """根据主机索引(0-191)确定其连接的TOR交换机ID(0-7)"""
    return host_idx // HOSTS_PER_TOR

# 生成作业文件
with open(f"{args.output_dir}/large_jobs_{args.tail_intensity}.txt", "w") as f:
    # 写入作业数量
    f.write(f"{NUM_JOBS}\n")
    
    for app_id in range(NUM_JOBS):
        # 为每个作业随机选择一个PS，确保PS的job数量不超过MAX_JOBS_PER_PS
        available_ps_indices = [i for i in range(NUM_HOSTS) if ps_job_count[i] < MAX_JOBS_PER_PS and not worker_allocation[i]]
        if not available_ps_indices:
            print(f"警告: 没有可用的PS主机了，只能生成{app_id}个作业")
            break
        
        ps_host_idx = random.choice(available_ps_indices)
        ps_job_count[ps_host_idx] += 1
        ps_node = get_host_real_index(ps_host_idx)  # 转换为实际节点ID
        
        # 随机选择worker数量（2-8个）
        worker_count = random.randint(2, 8)
        
        # 随机选择worker节点，确保不与已分配的worker节点冲突，且不选择PS节点
        available_workers = [i for i in range(NUM_HOSTS) if not worker_allocation[i] and i != ps_host_idx]
        if len(available_workers) < worker_count:
            print(f"警告: 可用worker节点不足，只能分配{len(available_workers)}个而不是{worker_count}个")
            worker_count = len(available_workers)
        
        selected_worker_indices = random.sample(available_workers, worker_count)
        
        # 标记已分配的worker节点
        for worker_idx in selected_worker_indices:
            worker_allocation[worker_idx] = True
        
        # 其他参数
        max_count = 100  # 每个worker发送的最大数据包数
        
        # # 使用tail intensity调整发送间隔（模拟拖尾）
        # base_interval = 0.00001  # 基准发送间隔
        # interval_variation = base_interval * tail_config['interval_variation']
        # interval = base_interval + random.uniform(-interval_variation, interval_variation)
        # interval = max(0.000001, interval)  # 确保间隔不会太小
        interval = 1.0
        # 确定聚合器节点 - 选择第一个worker所在的TOR交换机
        first_worker_idx = selected_worker_indices[0]
        worker_tor = get_host_tor(first_worker_idx)
        aggregator_node = worker_tor  # 使用worker连接的TOR交换机作为聚合器
        
        aggr_used = 450  # 聚合器使用的缓冲区大小
        model_size = 40  # 模型大小，单位MB
        start_time = 2.0 # 起始时间
        
        # 判断是否跨机架通信 - 检查worker是否连接到多个不同的TOR
        worker_tors = set(get_host_tor(widx) for widx in selected_worker_indices)
        cross_rack = 1 if len(worker_tors) > 1 or worker_tor != get_host_tor(ps_host_idx) else 0
        
        # 写入作业信息：appid, worker数量, 聚合器节点, PS节点, max_count, 间隔, aggr_used, 模型大小, 启动时间, cross_rack
        f.write(f"{app_id} {worker_count} {aggregator_node} {ps_node} {max_count} {interval} {aggr_used} {model_size} {start_time} {cross_rack}\n")
        
        # 写入每个worker节点
        for worker_idx in selected_worker_indices:
            real_worker_node = get_host_real_index(worker_idx)
            f.write(f"{real_worker_node} ")
        f.write("\n")

print(f"已生成{NUM_JOBS}个作业的配置文件，拖尾强度为{args.tail_intensity}")
print(f"作业文件保存在: {args.output_dir}/large_jobs_{args.tail_intensity}.txt")
print(f"配置文件保存在: {args.config_dir}/test[atp|a2tp|time].txt") 