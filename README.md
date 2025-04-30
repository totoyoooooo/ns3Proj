

### 说明

可以在命令中设置 `timewindow` 的大小，单位是秒，例如以下命令
```
./waf --run "testswitchml/start_test --configpath=./lzy_mix/config/testatp.txt --cmd_poolsize=450 --timewindow=0.0000001" 
```
就是设置时间窗口为 `1e-7s`

其中，`aggregator.cc` 中，`ForceFlush` 是定时器强制清空的函数，每次触发时会在终端输出形如 `ForceFlush called for appid: 0 key: 15615 after 1e-09 seconds` 的语句，可以借此判断是否触发了强制清空。

`aggregator.cc` 中，`SetOnoffTimeWindow` 函数的作用是切换为 `TimeWindow` 逻辑，目前代码是在 `2.000010s` 切换，如果不希望切换，可以注释掉 `StartApplication` 中的以下语句。

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