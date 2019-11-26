% test file
close all;clc; 
% strName = 'x .* sin(pi .* x) - y.* sin(pi .* y)';
% strName = 'x .* sin(pi .* x) - y.* sin(pi .* y)';
% funcStr = ['@(x,y)', strName];
% fh = str2func(funcStr);
% [x,y] = meshgrid(-1:0.05:1,-1:0.05:1);
% z = fh(x,y);
% mesh(x,y,z,'FaceColor',[0.5,0.5,0.5]);
% xlabel('x');ylabel('y');zlabel('f(x,y)');title(strName);

syms x y;
f = x .* sin(pi .* x) - y.* sin(pi .* y);
x0=[0.1;0.3];
e=10^(-20);
[k ender]=steepest(f,x0,e);
% plot3([x(1),ender], y,fx,'k*');


% 《神经网络模型及其matlab仿真程序设计》
%{
close all;clc;
% 一维数据,一列代表一个样本
[input,targ] = simplefit_dataset;
funcName = '2.*sin(x) + cos(y)';
% funcName = '-1 * x .* sin(2 * pi .* x) + y.* sin(2 * pi .* y) + 1'; 
funcName = 'sin(2 * pi .* x) + cos(2 * pi .* y) + 1'; 
funcStr = ['@(x, y)', funcName];
fh = str2func(funcStr);
[x,y] = meshgrid(-1:0.05:1,-1:0.05:1);
z = fh(x, y);

% 网络模型
targ = z(:)';
input = [x(:)';y(:)'];
% 隐藏层有10个神经元,
% feedforwardnet自动将输入数据分为训练\验证\测试子集
net = feedforwardnet(10, 'trainlm');
% 网络参数的设置
% net.divideParam.trainRatio = 0.8;
% net.divideParam.valRatio   = 0.1;
% net.divideParam.testRatio  = 0.1;

% net.trainParam.epochs=10000;    %训练次数设置
% net.trainParam.goal=1e-7;       %训练目标设置
% net.trainParam.lr=0.01;         %学习率设置
% net.trainParam.mc=0.9;          %动量因子的设置，默认为0.9
% net.trainParam.show=25;         %显示的间隔次数
% 增加隐藏层的数量

% 开始训练
net = train(net, input, targ);
% view(net)
outTest = net(input);
% perf = mse(y-t);
perf = perform(net, outTest, targ);
fprintf('perfor:%f\n', perf);

% 画图验证
mesh(x,y,z,'FaceColor','r');hold on;
mesh(x,y,reshape(outTest, size(x)),'FaceColor','b');
legend('target','test');
xlabel('x');ylabel('y');zlabel('f(x,y)');title(funcName);
%}