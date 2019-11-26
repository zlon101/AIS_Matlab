function [stegObj, nChanges, nzAC,k] = F5(cover, stegoPath, payload)
% -------------------------------------------------------------------------
% cover可以是路径或jobj结构体
% 根据嵌入效率计算出修改量, 嵌入效率考虑了缩减效应
% -------------------------------------------------------------------------

% 参数为文件路径
if( ischar(cover) )
    jobj = jpeg_read(cover);    % JPEG image structure
    % QDCT = jobj.coef_arrays{1};  % DCT plane    
    [QDCT,nChanges, nzAC,k] = modify(jobj.coef_arrays{1}, payload);
    jobj.coef_arrays{1} = QDCT;
    jobj.optimize_coding = 1;
    stegObj = jobj;
% 参数为jobj结构体
elseif( isstruct(cover) )
    [QDCT,nChanges, nzAC,k] = modify(cover.coef_arrays{1}, payload);
    stegObj = cover;
    stegObj.coef_arrays{1} = QDCT;
    stegObj.optimize_coding = 1;
% 参数为图像像素
else
    [QDCT, QTable, Cb, Cr] = getQDCT(cover, quality);
    [QDCT,nChanges, nzAC,k] = modify(QDCT, payload, isreduce);
    stegoPixel = qdct2Img(QDCT, QTable, Cb, Cr);
end
% changeRata = nChanges/numel(QDCT);
% fprintf('F5 changeRata:%d\n', changeRata);
jpeg_write(stegObj, stegoPath);
end

% 修改QDCT
function [QDCT,numChanged, numNzAC,k] = modify(QDCT, payload)
% payload:          嵌入率或者需要嵌入的数据量
seed = 99;rng('default');rng(seed);
if (payload > 0)    
	[numChanged, k, nEmbed] = getEfficiency(QDCT,payload);  % nEmbed:需要嵌入的数据量
    % numChanged：被修改系数的数量
    changeableInd = (QDCT~=0);
    changeableInd(1:8:end,1:8:end) = false;
    numNzAC = nnz(changeableInd);                     % 非零Ac系数的数量
    % 可修改的系数的位置
    changeableInd = find(changeableInd);              % indexes of the changeable coefficients
    %rand('state',seed);                              % initialize PRNG using given SEED    
    % 置乱
    changeableInd = changeableInd(randperm(numNzAC)); % create a pseudorandom walk over nonzero AC coefficients        
    to_be_changed = changeableInd(1:numChanged);      % coefficients to be changed    
    % 绝对值减一    
    QDCT(to_be_changed) = QDCT(to_be_changed) - sign(QDCT(to_be_changed));
end
end

function [numChanged,k,nEmbed,effci] = getEfficiency(Qdct,payload)
% numChanged:       被修改的数量
% payload:          嵌入率或者需要嵌入的数据量
% nEmbed：           需要嵌入的数据量
%{
计算嵌入效率k
switch(double(payload))
    case 0.6
        k=2;
    case 0.5
        k=2;
    case 0.4
        k=3;
    case 0.3
        k=3;
    case 0.2
        k=4;
    case 0.1
        k=5;
    otherwise
        k =1;
        fprintf('Unexpected payload:%f\n', payload);
        pause;
end
%}

% -------------------------------------------------------------------------
% 根据f5-java程序修改嵌入效率的计算公式
% 非零Ac系数, 用于计算实际嵌入量
numNzAc = nnz(Qdct~=0)-nnz(Qdct(1:8:end,1:8:end)~=0);
Qdct(1:8:end,1:8:end) = inf;
numZero = nnz(Qdct==0);         % 值为0的Ac个数
numOne = nnz(abs(Qdct)==1);     % 值为1的个数
numDc = nnz(Qdct==inf);         % DC系数的个数
large = numel(Qdct) - numZero - numOne - numDc;
capcity = large + 0.49*numOne;  % 嵌入容量
% 需要嵌入的数据量nEmbed
if(payload>5)
    nEmbed = payload;
else
    nEmbed = numNzAc*payload;       % payload: bpac-每个非零Ac系数携带多少信息
end

% 可嵌入数据量canEmbed
canEmbed = 0;
for k=7:-1:1                    % 注意k==1
    % k越小, canEmbed越大
    n = 2^k-1;
    canEmbed = (capcity/n)*k - mod((capcity/n)*k, n);
    if(canEmbed>=nEmbed)
        break;
    end
end

switch(k)
    case 1 
        effci=1.4;
    case 2 
        effci=1.6;
    case 3 
        effci=2.0;
    case 4 
        effci=2.3;
    case 5 
        effci=2.6;
    case 6 
        effci=3.2;
    case 7 
        effci=3.2;
end
%{
% 计算嵌入效率(结果与java计算的结果不同)
t = large - mod(large,n+1);
t = (t + numOne*1.5 - numOne/(n+1)) / (n+1);
effci=canEmbed/t + 0.1*mod( canEmbed*10/t, 10);
%}
% 被修改的数量
numChanged = round(nEmbed/effci);
end