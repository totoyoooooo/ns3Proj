#!/bin/bash

# 确保目录存在
mkdir -p results
mkdir -p lzy_mix/job
mkdir -p lzy_mix/config
mkdir -p lzy_mix/topology

# 协议列表
ALL_PROTOCOLS=("atp" "a2tp" "time")
# 拖尾强度列表
ALL_TAIL_INTENSITIES=("low" "medium" "high")
# 每个配置运行的轮数
NUM_ROUNDS=3

# 解析命令行参数
SPECIFIED_PROTOCOL=$1
SPECIFIED_TAIL=$2

# 如果指定了特定协议，则只运行该协议
if [ ! -z "$SPECIFIED_PROTOCOL" ]; then
    # 检查协议是否有效
    if [[ ! " ${ALL_PROTOCOLS[@]} " =~ " ${SPECIFIED_PROTOCOL} " ]]; then
        echo "错误: 无效的协议。有效的协议是: ${ALL_PROTOCOLS[@]}"
        exit 1
    fi
    PROTOCOLS=("$SPECIFIED_PROTOCOL")
else
    PROTOCOLS=("${ALL_PROTOCOLS[@]}")
fi

# 如果指定了特定拖尾强度，则只运行该拖尾强度
if [ ! -z "$SPECIFIED_TAIL" ]; then
    # 检查拖尾强度是否有效
    if [[ ! " ${ALL_TAIL_INTENSITIES[@]} " =~ " ${SPECIFIED_TAIL} " ]]; then
        echo "错误: 无效的拖尾强度。有效的拖尾强度是: ${ALL_TAIL_INTENSITIES[@]}"
        exit 1
    fi
    TAIL_INTENSITIES=("$SPECIFIED_TAIL")
else
    TAIL_INTENSITIES=("${ALL_TAIL_INTENSITIES[@]}")
fi

# 输出运行配置
echo "==== 大规模仿真运行配置 ===="
echo "协议: ${PROTOCOLS[@]}"
echo "拖尾强度: ${TAIL_INTENSITIES[@]}"
echo "每个配置运行轮数: $NUM_ROUNDS"
echo "============================"

# 创建一个用于存储结果的CSV文件
RESULTS_FILE="results/simulation_results.csv"
# 如果是完整运行则重新创建结果文件，否则只在不存在时创建
if [ ${#PROTOCOLS[@]} -eq ${#ALL_PROTOCOLS[@]} ] && [ ${#TAIL_INTENSITIES[@]} -eq ${#ALL_TAIL_INTENSITIES[@]} ]; then
    echo "Protocol,TailIntensity,Round,FinishTime,Throughput" > $RESULTS_FILE
else
    # 如果是部分运行且文件不存在，创建文件
    if [ ! -f "$RESULTS_FILE" ]; then
        echo "Protocol,TailIntensity,Round,FinishTime,Throughput" > $RESULTS_FILE
    fi
fi

# 运行模拟
for tail in "${TAIL_INTENSITIES[@]}"; do
    # 首先生成该拖尾强度的拓扑文件
    echo "生成拖尾强度为 ${tail} 的拓扑文件..."
    python3 generate_topology.py --tail_intensity "$tail"
    
    # 检查拓扑文件是否存在
    LARGE_TOPO_FILE="./lzy_mix/topology/topology_large_${tail}.txt"
    TEST_TOPO_FILE="./lzy_mix/topology/testtopo_${tail}.txt"
    
    if [ ! -f "$LARGE_TOPO_FILE" ]; then
        echo "错误: 大规模拓扑文件不存在: $LARGE_TOPO_FILE"
        exit 1
    fi
    
    if [ ! -f "$TEST_TOPO_FILE" ]; then
        echo "错误: 测试拓扑文件不存在: $TEST_TOPO_FILE"
        exit 1
    fi
    
    echo "拓扑文件已生成:"
    echo "  - 大规模拓扑: $LARGE_TOPO_FILE"
    echo "  - 测试拓扑: $TEST_TOPO_FILE"
    
    # 生成该拖尾强度的作业文件
    echo "生成拖尾强度为 ${tail} 的作业文件..."
    # python3 generate_jobs.py --tail_intensity "$tail" --jobs 10
    
    # 更新配置文件中的拓扑文件路径
    for protocol in "${PROTOCOLS[@]}"; do
        # 更新配置文件中的拓扑文件路径
        CONFIG_FILE="lzy_mix/config/test${protocol}.txt"
        
        # 检查配置文件是否存在
        if [ ! -f "$CONFIG_FILE" ]; then
            echo "错误: 配置文件不存在: $CONFIG_FILE"
            exit 1
        fi
        
        # 更新配置文件中的拓扑路径
        echo "更新配置文件 $CONFIG_FILE 中的拓扑路径..."
        sed -i "s|TOPO_FILE .*|TOPO_FILE ${LARGE_TOPO_FILE}|g" "$CONFIG_FILE"
        
        # 输出配置文件内容以进行验证
        echo "配置文件 $CONFIG_FILE 内容:"
        cat "$CONFIG_FILE"
        
        for round in $(seq 1 $NUM_ROUNDS); do
            echo "运行 ${protocol} 协议，拖尾强度: ${tail}，第 ${round} 轮..."
            
            # 运行模拟
            OUTPUT_FILE="results/${protocol}_${tail}_round${round}.txt"
            
            # 直接在命令行打印出完整命令，以便调试
            COMMAND="./waf --run \"testswitchml/start_test --configpath=${CONFIG_FILE} --cmd_poolsize=450 --timewindow=0.00000001 --tail=${tail} \""
            echo "执行命令: $COMMAND"
            
            # 执行命令
            ./waf --run "testswitchml/start_test --configpath=${CONFIG_FILE} --cmd_poolsize=450 --timewindow=0.00000001 --tail=${tail} " > "$OUTPUT_FILE" 2>&1
            
            # 检查命令执行结果
            if [ $? -ne 0 ]; then
                echo "错误: 模拟运行失败，检查 $OUTPUT_FILE 获取详细信息"
            else
                echo "模拟成功完成"
            fi
            
            # 提取完成时间（增强容错）
            # 从 "finish_time +X.XXXs total=Y" 格式中提取时间部分
            FINISH_TIME=$(grep "finish_time" "$OUTPUT_FILE" | tail -n 1 | awk -F'[ +]' '{if (NF >=2) print $3; else print "N/A"}' | sed 's/s$//')

            echo "提取到的完成时间: $FINISH_TIME"

            # 调试：打印吞吐量相关行
            echo "调试: 查看吞吐量相关行"
            grep -i "APP.*throughput" "$OUTPUT_FILE" | tail -5

            # 分析原始行格式确定列位置
            echo "分析第一行格式..."
            first_line=$(grep -i "APP.*throughput" "$OUTPUT_FILE" | head -1)
            echo "样本行: $first_line"
            # 提取 "aggregate throughput" 后面的数值列
            column_index=$(echo "$first_line" | awk '{for(i=1;i<=NF;i++) if($i=="throughput") print i+1}')
            echo "吞吐量值所在列: $column_index"

            # 提取稳定阶段(超过2秒)后的吞吐量值并计算平均值
            echo "尝试提取吞吐量值..."
            throughputs=$(grep -i "APP.*aggregate throughput" "$OUTPUT_FILE" | awk -v col="$column_index" '
            {
                # 提取时间（去掉+和s）
                time_str = $4;
                gsub(/\+|s/, "", time_str);
                time = time_str + 0;
                
                # 检查时间是否大于2秒（稳定阶段）且值有效
                if (time > 2.0 && $(col) != "" && $(col) ~ /^[0-9]+(\.[0-9]+)?$/) {
                    print $(col);
                }
            }')

            # 计算平均值（空数据时设为 N/A）
            THROUGHPUT_AVG=$(echo "$throughputs" | awk '
            {
                sum += $1;
                count++;
            } 
            END {
                if (count > 0) {
                    printf "%.4f", sum/count;
                } else {
                    print "N/A";  # 无有效数据时输出 N/A
                }
            }')

            # 输出结果
            echo "Stable Average Throughput: $THROUGHPUT_AVG Gbps"
            echo "提取到的吞吐量值列表: $throughputs"
            echo "有效吞吐量值数量: $(echo "$throughputs" | wc -l)"

            # 验证是否为有效数字
            if [[ ! "$THROUGHPUT_AVG" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                echo "警告: 未能从$OUTPUT_FILE中提取有效的吞吐量数据，将使用N/A"
                THROUGHPUT_AVG="N/A"
            fi

            # 记录结果
            echo "${protocol},${tail},${round},${FINISH_TIME:-N/A},${THROUGHPUT_AVG:-N/A}" >> "$RESULTS_FILE"
            # 等待几秒，确保系统资源释放
            sleep 2
        done
    done
done

# 分析结果 - 只在有结果数据的情况下进行
if [ -f "$RESULTS_FILE" ] && [ $(wc -l < "$RESULTS_FILE") -gt 1 ]; then
    echo "计算平均值和标准偏差..."
    python3 - << EOF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
import os

# 读取结果
try:

    df = pd.read_csv('$RESULTS_FILE')
    
    # 检查数据帧是否为空
    if df.empty:
        print("警告: 结果文件为空，无法进行分析")
        exit(0)
    
    # 检查是否有缺失值
    missing_finish = df['FinishTime'].isna().sum()
    missing_throughput = df['Throughput'].isna().sum()
    
    if missing_finish > 0 or missing_throughput > 0:
        print(f"警告: 结果中有 {missing_finish} 个缺失的完成时间和 {missing_throughput} 个缺失的吞吐量值")
    
    # 将N/A替换为NaN以便进行数值计算
    df.replace('N/A', np.nan, inplace=True)
    # 转换为数值类型
    df['FinishTime'] = pd.to_numeric(df['FinishTime'], errors='coerce')
    df['Throughput'] = pd.to_numeric(df['Throughput'], errors='coerce')
    
    # 计算每个协议和拖尾强度的平均值和标准偏差
    groups = df.groupby(['Protocol', 'TailIntensity'])
    summary = groups.agg({
        'FinishTime': ['mean', 'std'],
        'Throughput': ['mean', 'std']
    }).reset_index()
    
    # 保存汇总结果
    summary.to_csv('results/summary_results.csv')
    
    # 打印结果
    print("\n平均完成时间和吞吐量:")
    print(summary)
    
    # 确保结果目录存在
    os.makedirs('results', exist_ok=True)
    
    # 如果只有一种协议或一种拖尾强度，不需要绘制比较图表
    if len(df['Protocol'].unique()) == 1 and len(df['TailIntensity'].unique()) == 1:
        protocol = df['Protocol'].unique()[0]
        tail = df['TailIntensity'].unique()[0]
        avg_finish = df['FinishTime'].mean()
        avg_throughput = df['Throughput'].mean()
        print(f"\n协议 {protocol} 在 {tail} 拖尾强度下的性能:")
        print(f"平均完成时间: {avg_finish:.4f} 秒")
        print(f"平均吞吐量: {avg_throughput:.4f} Gbps")
        # 创建简单的轮次与性能的散点图
        plt.figure(figsize=(10, 5))
        plt.scatter(df['Round'], df['FinishTime'], label='完成时间(秒)', marker='o', color='blue')
        plt.title(f'协议 {protocol} 在 {tail} 拖尾强度下的完成时间')
        plt.xlabel('轮次')
        plt.ylabel('完成时间 (秒)')
        plt.grid(True)
        plt.savefig(f'results/{protocol}_{tail}_finish_time.png')
        
        plt.figure(figsize=(10, 5))
        plt.scatter(df['Round'], df['Throughput'], label='吞吐量(Gbps)', marker='^', color='green')
        plt.title(f'协议 {protocol} 在 {tail} 拖尾强度下的吞吐量')
        plt.xlabel('轮次')
        plt.ylabel('吞吐量 (Gbps)')
        plt.grid(True)
        plt.savefig(f'results/{protocol}_{tail}_throughput.png')
        exit(0)
    
    # 绘制图表
    
    # 颜色和标记的映射
    protocol_colors = {'atp': 'blue', 'a2tp': 'red', 'time': 'green'}
    tail_markers = {'low': 'o', 'medium': 's', 'high': '^'}
    
    # 设置图表样式
    plt.style.use('ggplot')
    plt.rcParams['figure.figsize'] = (12, 7)
    plt.rcParams['font.size'] = 12
    
    # 如果有多种协议，绘制协议比较图
    if len(df['Protocol'].unique()) > 1:
        # 1. 每种协议的平均完成时间对比图
        plt.figure()
        for protocol in df['Protocol'].unique():
            data = summary[summary['Protocol'] == protocol]
            if len(data) > 0:  # 确保有数据
                plt.errorbar(
                    data['TailIntensity'], 
                    data[('FinishTime', 'mean')], 
                    yerr=data[('FinishTime', 'std')], 
                    label=protocol.upper(),
                    marker='^',
                    markersize=8,
                    capsize=5,
                    color=protocol_colors.get(protocol, 'black'),
                    linewidth=2
                )
        
        plt.title('各协议在不同拖尾强度下的平均完成时间', fontsize=16)
        plt.xlabel('拖尾强度', fontsize=14)
        plt.ylabel('完成时间 (秒)', fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('results/protocols_finish_time_comparison.png')
        
        # 2. 每种协议的平均吞吐量对比图
        plt.figure()
        for protocol in df['Protocol'].unique():
            data = summary[summary['Protocol'] == protocol]
            if len(data) > 0:  # 确保有数据
                plt.errorbar(
                    data['TailIntensity'], 
                    data[('Throughput', 'mean')], 
                    yerr=data[('Throughput', 'std')], 
                    label=protocol.upper(),
                    marker='o',
                    markersize=8,
                    capsize=5,
                    color=protocol_colors.get(protocol, 'black'),
                    linewidth=2
                )
        
        plt.title('各协议在不同拖尾强度下的平均吞吐量', fontsize=16)
        plt.xlabel('拖尾强度', fontsize=14)
        plt.ylabel('吞吐量 (Gbps)', fontsize=14)
        plt.legend(fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('results/protocols_throughput_comparison.png')
    
    # 如果有多种拖尾强度，绘制拖尾强度比较图
    if len(df['TailIntensity'].unique()) > 1:
        # 3. 各拖尾强度下不同协议的完成时间对比
        fig, axes = plt.subplots(1, len(df['TailIntensity'].unique()), figsize=(18, 6))
        tail_intensities = sorted(df['TailIntensity'].unique())
        
        # 处理只有一种拖尾强度的情况（确保axes是列表）
        if len(tail_intensities) == 1:
            axes = [axes]
        
        for i, tail in enumerate(tail_intensities):
            tail_data = df[df['TailIntensity'] == tail]
            for protocol in df['Protocol'].unique():
                protocol_data = tail_data[tail_data['Protocol'] == protocol]
                if len(protocol_data) > 0:  # 确保有数据
                    axes[i].scatter(
                        protocol_data['Round'], 
                        protocol_data['FinishTime'], 
                        label=protocol.upper(),
                        marker=tail_markers.get(tail, 'o'),
                        s=100,
                        color=protocol_colors.get(protocol, 'black'),
                        alpha=0.7
                    )
            
            axes[i].set_title(f'拖尾强度: {tail}', fontsize=14)
            axes[i].set_xlabel('轮次', fontsize=12)
            axes[i].set_ylabel('完成时间 (秒)', fontsize=12)
            axes[i].legend()
            axes[i].grid(True)
        
        plt.tight_layout()
        plt.savefig('results/tail_intensity_finish_time_comparison.png')
        
        # 4. 各拖尾强度下不同协议的吞吐量对比
        fig, axes = plt.subplots(1, len(df['TailIntensity'].unique()), figsize=(18, 6))
        
        # 处理只有一种拖尾强度的情况
        if len(tail_intensities) == 1:
            axes = [axes]
        
        for i, tail in enumerate(tail_intensities):
            tail_data = df[df['TailIntensity'] == tail]
            for protocol in df['Protocol'].unique():
                protocol_data = tail_data[tail_data['Protocol'] == protocol]
                if len(protocol_data) > 0:  # 确保有数据
                    axes[i].scatter(
                        protocol_data['Round'], 
                        protocol_data['Throughput'], 
                        label=protocol.upper(),
                        marker=tail_markers.get(tail, 'o'),
                        s=100,
                        color=protocol_colors.get(protocol, 'black'),
                        alpha=0.7
                    )
            
            axes[i].set_title(f'拖尾强度: {tail}', fontsize=14)
            axes[i].set_xlabel('轮次', fontsize=12)
            axes[i].set_ylabel('吞吐量 (Gbps)', fontsize=12)
            axes[i].legend()
            axes[i].grid(True)
        
        plt.tight_layout()
        plt.savefig('results/tail_intensity_throughput_comparison.png')
    
    # 5. 每种协议随拖尾强度变化的性能降低百分比（只在有低拖尾数据的情况下计算）
    if 'low' in df['TailIntensity'].unique() and len(df['TailIntensity'].unique()) > 1:
        plt.figure(figsize=(12, 7))
        
        # 计算每种协议在低拖尾强度下的基准性能
        baseline = df[df['TailIntensity'] == 'low'].groupby('Protocol')['Throughput'].mean()
        
        # 计算中、高拖尾强度下相对于基准的性能变化百分比
        for tail in ['medium', 'high']:
            if tail in df['TailIntensity'].unique():
                performance = df[df['TailIntensity'] == tail].groupby('Protocol')['Throughput'].mean()
                percentage = ((performance - baseline) / baseline) * 100
                
                # 绘制柱状图
                plt.bar(
                    [p + ' (' + tail + ')' for p in percentage.index], 
                    percentage.values,
                    alpha=0.7,
                    label=f'{tail} 拖尾'
                )
        
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        plt.title('各协议在不同拖尾强度下的吞吐量变化百分比（相对于低拖尾）', fontsize=16)
        plt.xlabel('协议 (拖尾强度)', fontsize=14)
        plt.ylabel('吞吐量变化百分比 (%)', fontsize=14)
        plt.grid(True, axis='y')
        plt.legend()
        plt.tight_layout()
        plt.savefig('results/performance_degradation_percentage.png')
    
    print("\n分析完成! 图表已保存到 results 目录。")
except Exception as e:
    print(f"分析时发生错误: {str(e)}")
    import traceback
    traceback.print_exc()
EOF
fi

# 根据运行配置输出不同的完成消息
if [ ${#PROTOCOLS[@]} -eq 1 ] && [ ${#TAIL_INTENSITIES[@]} -eq 1 ]; then
    echo "仿真完成! 协议 ${PROTOCOLS[0]} 在 ${TAIL_INTENSITIES[0]} 拖尾强度下的结果已保存在 results 目录中。"
else
    echo "仿真完成! 所有结果保存在 results 目录中。"
fi 