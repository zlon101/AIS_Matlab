function [bestFits,bestAbs,meanFits,Memory] = CSA(srcPath)
% 克隆选择算法, 对图像锐化参数进行优化
% 
%%
Root = 'E:\astego\Images\Experis\';
% name = split(srcPath, '\');  name = name(end);
name = 'xx.pgm';
srcStegoPath= [Root,'stegos\', name];
sharpedPath = [Root,'sharpeds\',name];
sharpedStegoPath = [Root,'sharpedStegos\',name];
payLoad = single(0.4);
embedParas = struct('srcPath',srcPath,'sharpedPath',sharpedPath,...
    'sharpedStegoPath',sharpedStegoPath,'payLoad',payLoad);
% srcData = single(imread(srcPath));
% srcStegoData = HUGO_like(uint8(srcData), payLoad);
% imwrite(uint8(srcStegoData),srcStegoPath, 'pgm');

Precision = 0.001;
Vmin = 0;  Vmax = 1.4;
L = log2( ((Vmax-Vmin)/Precision) + 1);
L = ceil(L);  % 编码长度
NumParas = 1;  % 参数个数
NumTotal = single(15);  % 抗体个数
Iters = 8;  % 迭代次数
Memory = containers.Map('KeyType','char','ValueType','double');
% output
bestFits = zeros(Iters, 1);
bestAbs = cell(Iters, 1);
meanFits = zeros(Iters,1);
% 初始化
MultRata = 0.3;
PClone = 0.3;  NumCloned = round(NumTotal * PClone);
PMuMin = 0.02;  PMuMax = 0.1;  PMu = PMuMin;
PNewMin = 0.1; PNewMax = 0.3; PNew= NumTotal*PNewMin;
T = 6;
genes = initAb(NumTotal, NumParas*L);
% 0值编码
genes(1,:) = zeros(1,NumParas*L);
% ------------------测试Castro---------------------------------
%{
f = '1 * x .* sin(4 * pi .* x) - 1 * y.* sin(4 * pi .* y + pi) + 1';
[x,y] = meshgrid(-1:0.05:1,-1:0.05:1); vxp = x; vyp = y;
vzp = eval(f);  % 函数值
NumParas=2;
Abs = decodeAbs(genes, NumParas, Vmin, Vmax);
Abs = cell2mat(Abs);
x = Abs(:,1);
y = Abs(:,2);
fit = eval(f);
imprime(1,vxp,vyp,vzp,x,y,fit,1,1); title('Initial Population');
% ------------------测试end---------------------------------
%}

%% 开始迭代
for i=1:Iters
%% 计算适应度
Abs = decodeAbs(genes, NumParas,Vmin,Vmax);  % N*NumParas cell
[fits, Memory] = calcuFit(Abs, embedParas, Memory);
[fits, sortInd]= sort(fits, 'ascend');  % descend:降序, 要求优秀的排在前面
Abs= Abs(sortInd, :);
genes= genes(sortInd, :);

% 消亡
if(size(genes,1) > NumTotal)
    genes(NumTotal+1:end,:) = []; Abs(NumTotal+1:end,:)=[];
    fits(NumTotal+1:end,:)=[];
end
% 日志
bestFits(i) = fits(1);
bestAbs(i) = Abs(1);  % Abs: cell
meanFits(i) = mean(fits);
fprintf('\nIter:%d - best fit: %5.3f\n', i,fits(1));

if(fits(1) <=bestFits(end))
    countBreak = countBreak+1;
else
    countBreak = 1;
    PMu = PMuMin;
    PNew = PNewMin;
end
if(countBreak > T)
    break;
elseif(countBreak > 0.5*T)
    PMu = PMuMax;
    PNew = PNewMax;
end
%% 克隆
[tmpGenes, pcs] = reprod(genes, NumCloned, MultRata);
% 变异
M = rand(size(tmpGenes)) <= PMu;  % M=1, 0,1翻转, 否则不变
tmpGenes = tmpGenes - 2 .* (tmpGenes.*M) + M;
% 维持现有最优Ab
tmpGenes(pcs,:) = genes(1:NumCloned, :);
genes = [tmpGenes; initAb(NumTotal*PNew, NumParas*L)];
% for-end
end
end