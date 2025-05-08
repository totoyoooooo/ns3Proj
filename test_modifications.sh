#!/bin/bash

# 测试脚本 - 验证修改是否正确

# 1. 测试生成拓扑文件
echo "==== 测试生成拓扑文件 ===="
for intensity in low medium high; do
    echo "生成 $intensity 拖尾强度的拓扑..."
    ./generate_topology.py --tail_intensity $intensity
    
    # 检查文件是否存在
    if [ -f "./lzy_mix/topology/topology_large_${intensity}.txt" ] && [ -f "./lzy_mix/topology/testtopo_${intensity}.txt" ]; then
        echo "✓ 成功: 拓扑文件已生成"
        # 检查文件大小是否合理
        size1=$(stat -c%s "./lzy_mix/topology/topology_large_${intensity}.txt")
        size2=$(stat -c%s "./lzy_mix/topology/testtopo_${intensity}.txt")
        if [ $size1 -gt 1000 ] && [ $size2 -gt 500 ]; then
            echo "✓ 成功: 文件大小合理 ($size1 字节, $size2 字节)"
        else
            echo "✗ 失败: 文件大小不合理 ($size1 字节, $size2 字节)"
        fi
    else
        echo "✗ 失败: 拓扑文件未生成"
    fi
    echo
done

# 2. 测试生成作业配置
echo "==== 测试生成作业配置 ===="
for intensity in low medium high; do
    echo "生成 $intensity 拖尾强度的作业配置..."
    ./generate_jobs.py --tail_intensity $intensity --jobs 5
    
    # 检查作业文件是否存在
    if [ -f "./lzy_mix/job/large_jobs_${intensity}.txt" ]; then
        echo "✓ 成功: 作业文件已生成"
        # 检查文件大小是否合理
        size=$(stat -c%s "./lzy_mix/job/large_jobs_${intensity}.txt")
        if [ $size -gt 100 ]; then
            echo "✓ 成功: 文件大小合理 ($size 字节)"
        else
            echo "✗ 失败: 文件大小不合理 ($size 字节)"
        fi
    else
        echo "✗ 失败: 作业文件未生成"
    fi
    
    # 检查配置文件是否存在
    for protocol in atp a2tp time; do
        if [ -f "./lzy_mix/config/test${protocol}.txt" ]; then
            echo "✓ 成功: ${protocol} 配置文件已生成"
            
            # 检查配置文件中的拓扑路径是否正确
            if grep -q "TOPO_FILE ./lzy_mix/topology/topology_large_${intensity}.txt" "./lzy_mix/config/test${protocol}.txt"; then
                echo "✓ 成功: ${protocol} 配置文件中包含正确的拓扑路径"
            else
                echo "✗ 失败: ${protocol} 配置文件中拓扑路径不正确"
            fi
        else
            echo "✗ 失败: ${protocol} 配置文件未生成"
        fi
    done
    echo
done

echo "所有测试完成!"
 