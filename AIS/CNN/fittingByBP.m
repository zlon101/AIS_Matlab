% fittingByBP.m
%% 用神经网络(BP)进行数据拟合

% 1.newff：神经网络参数设置函数
% net=newff(P,T,S,TF,BTF,PF,IPF,OPF,DDF)
    % P:输入数据矩阵。 T：输出数据矩阵。 S:隐含层节点数。 TF:节点传递函数
    % BTF：训练函数，包括梯度下降BP算法训练函数traingd, 
        % 动量反传的梯度下降BP算法traingdm,动态自适应学习率的梯度下降BP算法训练函数traingda,
        % 动量反传和动态自适应学习率的梯度下降BP算法训练函数traingdx, Levenberg_Marquardt的BP算法训练函数trainlm。
    % BLF:网络学习函数
    % PF：性能分析函数，包括均值绝对误差mae，均方差mse。
    % IPF：输入处理函数，OPF：输出处理函数，DDF：验证数据划分函数
% feedforwardnet()

% 2.train:BP神经网络训练函数
% [net,tr]=train(NET, X,T,Pi, Ai)
    % NET:待训练网络。X：输入数据矩阵。T：输出数据矩阵。
    % Pi:初始化输入层条件。Ai:初始化输出层条件。
    % net:训练好的网络。tr:训练过程记录
    
% 3.sim:BP神经网络预测函数
% y=sim(net,x)
    % x:输入数据。y:网络预测数据。
%% ================================================================
close all;clc; 
strName = 'x .* sin(pi .* x) - y.* sin(pi .* y)';
funcStr = ['@(x,y)', strName];
fh = str2func(funcStr);
[x,y] = meshgrid(-1:0.05:1,-1:0.05:1);
z = fh(x,y);
mesh(x,y,z,'FaceColor',[0.5,0.5,0.5]);
xlabel('x');ylabel('y');zlabel('f(x,y)');title(strName);

%{
vy = y(:)';
dataInput = [x(:)';y(:)']; 
dataTarg = z(:)';
% 归一化映射(行向量)
dataInput = mapminmax(dataInput);
vy = mapminmax(vy);
dataTarg = mapminmax(dataTarg);
% 总的样本个数
numSamp = size(dataInput,2);
% 默认配置
paramConfi.Q = numSamp;
paramConfi.trainRatio = 0.7;
paramConfi.valRatio   = 0.15;
paramConfi.testRatio  = 0.15;
% 划分数据集(训练\验证\测试)
[paramConfi.trainInd, paramConfi.verInd, paramConfi.testInd] = ...
    dividerand(paramConfi.Q, paramConfi.trainRatio, paramConfi.valRatio, paramConfi.testRatio);
in.train = dataInput(:, paramConfi.trainInd);
in.ver = dataInput(:, paramConfi.verInd);
in.test = dataInput(:, paramConfi.testInd);
targ.train = dataTarg(:, paramConfi.trainInd);
targ.ver = dataTarg(:, paramConfi.verInd);
targ.test = dataTarg(:, paramConfi.testInd);

% 创建网络
TF1='tansig';TF2='purelin';
net=newff(in.train, targ.train, 1000, {TF1 TF2},'traingdm');
% 网络参数的设置
net.trainParam.epochs=10000;    %训练次数设置
net.trainParam.goal=1e-7;       %训练目标设置
net.trainParam.lr=0.01;         %学习率设置
net.trainParam.mc=0.9;          %动量因子的设置，默认为0.9
net.trainParam.show=25;         %显示的间隔次数
% 指定训练参数
net.trainFcn = 'traingd';       %梯度下降算法
% net.trainFcn = 'traingdm';      %动量梯度下降算法
% net.trainFcn = 'trainlm';
% 开始训练
[net,tr]=train(net, in.train, targ.train);

% 计算仿真，其一般用sim函数
[outVer,trainPerf]=sim(net, in.ver);
mse(outVer-targ.ver);
% 验证的数据，经BP得到的结果
% [normvalidateoutput,validatePerf]=sim(net,valsample.p,[],[],valsample.t);
% 测试数据，经BP得到的结果
% [normtestoutput,testPerf]=sim(net,testsample.p,[],[],testsample.t);

% 将所得的结果进行反归一化，得到其拟合的数据
% trainoutput=mapminmax('reverse',normtrainoutput,ts);
% validateoutput=mapminmax('reverse',normvalidateoutput,ts);
% testoutput=mapminmax('reverse',normtestoutput,ts);

% 正常输入的数据的反归一化的处理，得到其正式值
trainvalue=mapminmax('reverse',trainsample.t,ts);   % 正常的验证数据
validatevalue=mapminmax('reverse',valsample.t,ts);  % 正常的验证的数据
testvalue=mapminmax('reverse',testsample.t,ts);     % 正常的测试数据

% 测试
pnew=[313,256,239]';
pnewn=mapminmax(pnew);
anewn=sim(net,pnewn);
anew=mapminmax('reverse',anewn,ts);
% 绝对误差的计算
errors=trainvalue-trainoutput;
% plotregression拟合图
figure,plotregression(trainvalue,trainoutput)
% 误差图
figure,plot(1:length(errors),errors,'-b')
title('误差变化图')
% 误差值的正态性的检验
figure,hist(errors);%频数直方图
figure,normplot(errors);%Q-Q图
[muhat,sigmahat,muci,sigmaci]=normfit(errors); %参数估计 均值,方差,均值的0.95置信区间,方差的0.95置信区间
[h1,sig,ci]= ttest(errors,muhat);%假设检验
figure, ploterrcorr(errors);%绘制误差的自相关图
figure, parcorr(errors);%绘制偏相关图
%}
%{
k=rand(1,2000);
[m,n]=sort(k);
input_train=input(n(1:1900),:)';
output_train=output(n(1:1900),:)';
input_test=input(n(1901:2000),:)';
output_test=output(n(1901:2000),:);
[inputn,inputps]=mapminmax(input_train);
[outputn,outputps]=mapminmax(output_train);
%BP神经网络构建
net=newff(inputn,outputn,5);
%网络参数配置（迭代次数，学习率，目标）
net.trainParam.epochs=100;
net.trainParam.lr=0.1;
net.trainParam.goal=0.00004;

%神经网络训练
net=train(net,inputn,outputn);

%预测输入数据归一化
inputn_test=mapminmax('apply',input_test,inputps);
%神经网络预测输出
an=sim(net,inputn_test);
%输出反归一化
BPoutput=mapminmax('reverse',an,outputps);

%网络预测结果图形
figure(1)
plot(BPoutput,':og')
hold on
plot(output_test,'-*')
legend('预测输出','期望输出')
title('BP网络预测输出','fontsize',12)
ylabel('函数输出','fontsize',12)
xlabel('样本','fontsize',12)

%网络预测误差图形
figure(2)
plot(error,'-*')
title('BP预测误差','fontsize',12)
ylabel('误差',fontsize,12)
ylabel('样本',fontsize,12)
%}