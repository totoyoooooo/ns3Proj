#!/bin/bash

# 测试单一配置的脚本，用于调试
# 使用: ./test_single_config.sh [protocol] [tail_intensity]
# 例如: ./test_single_config.sh atp medium

# 默认值
PROTOCOL=${1:-"time"}
TAIL_INTENSITY=${2:-"low"}

# 检查参数
if [[ ! "$PROTOCOL" =~ ^(atp|a2tp|time)$ ]]; then
    echo "错误: 协议必须是 atp、a2tp 或 time"
    exit 1
fi

if [[ ! "$TAIL_INTENSITY" =~ ^(low|medium|high)$ ]]; then
    echo "错误: 拖尾强度必须是 low、medium 或 high"
    exit 1
fi

echo "===== 测试单一配置 ====="
echo "协议: $PROTOCOL"
echo "拖尾强度: $TAIL_INTENSITY"

# 确保目录存在
mkdir -p results
mkdir -p lzy_mix/job
mkdir -p lzy_mix/config
mkdir -p lzy_mix/topology

# 生成拓扑文件
echo "生成拖尾强度为 ${TAIL_INTENSITY} 的拓扑文件..."
python3 generate_topology.py --tail_intensity "$TAIL_INTENSITY"

# 检查拓扑文件是否存在
LARGE_TOPO_FILE="./lzy_mix/topology_large_${TAIL_INTENSITY}.txt"
TEST_TOPO_FILE="./lzy_mix/topology/testtopo_${TAIL_INTENSITY}.txt"

if [ ! -f "$LARGE_TOPO_FILE" ]; then
    echo "错误: 大规模拓扑文件不存在: $LARGE_TOPO_FILE"
    exit 1
fi

if [ ! -f "$TEST_TOPO_FILE" ]; then
    echo "错误: 测试拓扑文件不存在: $TEST_TOPO_FILE"
    exit 1
fi

# 显示拓扑文件的内容前几行，用于调试
echo "拓扑文件内容预览:"
echo "大规模拓扑文件 ($LARGE_TOPO_FILE):"
head -n 5 "$LARGE_TOPO_FILE"
echo "测试拓扑文件 ($TEST_TOPO_FILE):"
head -n 5 "$TEST_TOPO_FILE"

echo "拓扑文件已生成:"
echo "  - 大规模拓扑: $LARGE_TOPO_FILE"
echo "  - 测试拓扑: $TEST_TOPO_FILE"

# 生成作业文件
echo "生成拖尾强度为 ${TAIL_INTENSITY} 的作业文件..."
python3 generate_jobs.py --tail_intensity "$TAIL_INTENSITY" --jobs 20

# 更新配置文件中的拓扑文件路径
CONFIG_FILE="lzy_mix/config/test${PROTOCOL}.txt"

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

# 添加调试输出到日志文件
OUTPUT_FILE="results/${PROTOCOL}_${TAIL_INTENSITY}_debug.txt"

# 显示完整命令，方便调试
COMMAND="./waf --run \"testswitchml/start_test --configpath=${CONFIG_FILE} --cmd_poolsize=450 --timewindow=0.00000001 --tail=${TAIL_INTENSITY} --topology=${TEST_TOPO_FILE}\""
echo "执行命令: $COMMAND"

# 添加与拓扑文件有关的调试输出
echo "拓扑文件信息:"
echo "- 文件路径: ${TEST_TOPO_FILE}"
echo "- 文件大小: $(du -h ${TEST_TOPO_FILE} | cut -f1) ($(wc -l < ${TEST_TOPO_FILE}) 行)"
echo "- 文件权限: $(ls -l ${TEST_TOPO_FILE})"

# 运行仿真，捕获详细输出
echo "正在运行模拟，输出将保存到 $OUTPUT_FILE..."
./waf --run "testswitchml/start_test --configpath=${CONFIG_FILE} --cmd_poolsize=450 --timewindow=0.00000001 --tail=${TAIL_INTENSITY} --topology=${TEST_TOPO_FILE}" > "$OUTPUT_FILE" 2>&1

# 检查命令执行结果
if [ $? -ne 0 ]; then
    echo "错误: 模拟运行失败，检查 $OUTPUT_FILE 获取详细信息"
    # 显示日志文件的最后20行
    echo "日志文件的最后20行:"
    tail -n 20 "$OUTPUT_FILE"
else
    echo "模拟成功完成！"
    
    # 从输出文件中提取完成时间和吞吐量
    FINISH_TIME=$(grep "finish_time" "$OUTPUT_FILE" | awk '{print $2}')
    THROUGHPUT=$(grep "aggreate throughput" "$OUTPUT_FILE" | awk '{print $4}')
    
    # 如果没有找到吞吐量数据，则尝试另一种方式提取
    if [ -z "$THROUGHPUT" ]; then
        THROUGHPUT=$(grep "ans:" "$OUTPUT_FILE" | awk '{print $6}')
    fi
    
    # 显示结果
    echo "完成时间: ${FINISH_TIME:-未找到}"
    echo "吞吐量: ${THROUGHPUT:-未找到} Mbps"
    
    # 检查是否生成了特定协议的输出文件
    PROTOCOL_OUTPUT="PS_${PROTOCOL}.txt"
    if [ -f "$PROTOCOL_OUTPUT" ]; then
        echo "协议输出文件已生成: $PROTOCOL_OUTPUT"
        # 显示文件的前几行
        echo "文件内容预览:"
        head -n 5 "$PROTOCOL_OUTPUT"
    else
        echo "警告: 未找到协议输出文件 $PROTOCOL_OUTPUT"
    fi
fi

echo "测试完成！结果保存在 $OUTPUT_FILE" 