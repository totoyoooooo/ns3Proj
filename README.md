# ns3proj

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