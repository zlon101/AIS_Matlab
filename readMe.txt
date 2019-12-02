# 编译说明

M文件编译为MEX：coder
M文件编译为EXE: 
  1. mcc -m **.mcc
  2. App-->App Compiler
  
编译环境:
  minGW-w64
  SDK 10
  
运行环境:
  Matlab Runtime
  
# 环境部署
  mbuild -setup 提示安装minGW 和 SDK,
  根据链接下载 mingw.mlpkginstall文件,并在MATLAB中打开,安装MinGW-w64,
  根据链接下载winsdksetup.exe 或 ios文件,然后安装SDK
  
  安装Matlab Runtime: D:\MATLAB2017B\toolbox\compiler\deploy\win64\MCRInstaller.exe
  
# 打包部署
  for_redistribution: 安装包
  for_testing: 与MCC命令编译结果相同
  
  