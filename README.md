
### 运行大规模仿真


运行所有协议和拖尾强度的组合
```
./run_large_scale_sim.sh
```
仅运行特定协议和拖尾强度的组合
```
./run_large_scale_sim.sh atp low
./run_large_scale_sim.sh a2tp medium
./run_large_scale_sim.sh time high
```

所有的结果将保存在 results 目录中，包括原始输出、CSV 格式的统计数据以及可视化图表，便于比较不同协议在不同拖尾强度下的性能。

### 文件说明

目前设置5个任务，如果需要修改任务数量，请在`test_modifications.sh`中 修改以下语句
```
./generate_jobs.py --tail_intensity $intensity --jobs 5
```
然后重新运行`test_modifications.sh`


/job 中的这些文件是在不同拖尾强度(tail intensity)下生成的作业配置文件：
1. large_jobs_low.txt - 低拖尾强度下的作业配置
2. large_jobs_medium.txt - 中等拖尾强度下的作业配置
3. large_jobs_high.txt - 高拖尾强度下的作业配置
每个文件的格式如下：
- 第一行是作业数量（例如，large_jobs_low.txt中有20个作业，而large_jobs_medium.txt和large_jobs_high.txt各有5个作业）
- 后续每个作业的格式为： 
- 首行：app_id worker数量 聚合器节点 PS节点 max_count 间隔 aggr_used 模型大小 启动时间 cross_rack 
- 随后几行：参与此作业的worker节点ID列表
例如，large_jobs_high.txt第一个作业：
```
0 5 1 148 100 1.0 450 40 2 1
52 82 203 69 130
```
这表示：
- 作业ID为0
- 有5个worker节点
- 聚合器节点是交换机1
- PS节点是主机148
- 每个worker发送100个数据包
- 发送间隔约为1.0
- 聚合器使用450单位的缓冲区
- 模型大小为40MB
- 启动时间为2.0秒
- cross_rack=1（表示跨机架通信）

PS节点识别作为PS的主机在每个作业配置中明确指定。以large_jobs_high.txt为例：
- 作业0的PS节点是主机148
- 作业1的PS节点是主机174
- 作业2的PS节点是主机204
- 作业3的PS节点是主机80
- 作业4的PS节点是主机20

在整个网络拓扑中，主机节点的ID从16开始（前16个ID是交换机），所以所有ID大于等于16的节点都是主机。

每个主机连接到特定的TOR交换机，TOR交换机再连接到所有核心交换机。PS节点是作为参数传输模型的服务器主机，负责接收和处理所有worker节点发送的数据。对于每个作业，都有一个指定的主机作为PS节点，与多个worker节点协作完成分布式训练任务。

--------------------

### 说明

可以在命令中设置 `timewindow` 的大小，单位是秒，例如以下命令
```
./waf --run "testswitchml/start_test --configpath=./lzy_mix/config/testatp.txt --cmd_poolsize=450 --timewindow=0.0000001" 
```
就是设置时间窗口为 1e-7s

其中，aggregator.cc 中，ForceFlush是定时器强制清空的函数，每次触发时会在终端输出形如 `ForceFlush called for appid: 0 key: 15615 after 1e-09 seconds` 的语句，可以借此判断是否触发了强制清空。

aggregator.cc 中，SetOnoffTimeWindow 函数的作用是切换为 TimeWindow 逻辑，目前代码是在 2.000010s 切换，如果不希望切换，可以注释掉 StartApplication 中的以下语句。

```
Simulator::Schedule(Seconds(1.000010), &UdpAggregator::SetOnoffTimeWindow, this);
```
如果需要改变切换的时间，修改上述语句中的 `1.000010`，希望在开始发送后的 `X s` 时刻切换，就改为 `1 + X`



`*_aggr_rate.log` 存储各协议聚合率，比如 `ATP_aggr_rate.log` 就是 ATP 的聚合情况

`turnover_rate_*` 存储各协议下，梯度包占用聚合器的情况，比如 `turnover_rate_ATP`

这些文件都在 `NS3Proj` 目录下



### 编译构建
设置权限为可执行
```
chmod +x ./waf
```
configure
```
CXXFLAGS=-Wno-error CC=g++-10 GCC=g++-10 CXX=g++-10 ./waf configure
```
build
```
sudo ./waf build
```

### 测试命令
测 `timewindow`
```
./waf --run "testswitchml/start_test --configpath=./lzy_mix/config/testtime.txt --cmd_poolsize=250"
```

测 `atp`
```
./waf --run "testswitchml/start_test --configpath=./lzy_mix/config/testatp.txt --cmd_poolsize=250"
```

测 `a2tp`
```
./waf --run "testswitchml/start_test --configpath=./lzy_mix/config/testa2tp.txt --cmd_poolsize=250"
```

测 `switchml`
```
./waf --run "testswitchml/start_test --configpath=./lzy_mix/config/testswitchml.txt --cmd_poolsize=250"
```

### 切换为写入
在 `udp-aggregator.cc` 中，注释掉 `StartApplication` 中的
```
LoadCachedSamples();
```
取消 `StartApplication` 中的以下注释 
```
// m_timeLogFile.open(m_timeDataFile, std::ios::trunc);
```
取消 `HandleRead` 中的以下注释 
```
// LogArrivalTime(app_id, recv_key);  
```

### 切换为读取 sample.txt
在 `udp-aggregator.cc` 中，取消 `StartApplication` 中的以下注释
```
// LoadCachedSamples();
```
注释掉 `StartApplication` 中的以下内容
```
m_timeLogFile.open(m_timeDataFile, std::ios::trunc);
```
注释掉 `HandleRead` 中的以下内容 
```
LogArrivalTime(app_id, recv_key);  
```

### 流量监测图
```
chmod +x testswitchml/run_model_comparison.sh
./testswitchml/run_model_comparison.sh
```
编译成功后，会弹出以下选择
```
Which model would you like to run? (1 for ResNet, 2 for VGG, 3 for both)
```
之后输入对应的数字来运行对应代码，得到对应图像

生成的图像在 `plots` 文件夹下