function QDCT = changeLFQDCT(QDCT, indSubBlockLF, changeRate, changeRateZero, debug)
% 修改低频QDCT系数
% indSubBlockLF:一个子块的低频系数索引
% changeRate:      值为1的系数修改概率
% changeRateZero   值为0的系数修改概率
% indSubBlockLF = [8,8; 8,7; 7,8; 8,6;7,7;6,8;8,5;7,6;6,7;5,8];

% 所选的低频系数的个数
nLF = size(indSubBlockLF, 1);
rQDCT = size(QDCT,1);
cQDCT = size(QDCT,2);
nBlock = floor(size(QDCT,1)/8) * floor(size(QDCT,2)/8);
% 子块二维索引
nRR = floor(size(QDCT,1)/8);
nCC = floor(size(QDCT,2)/8);
% 第一列表示第一维的索引
indBlock = zeros(nBlock,2);
k=1;
for i=1:nRR
    indBlock(k:k+nCC-1, 1) = i;
    indBlock(k:k+nCC-1, 2) = 1:nCC;
    k = k+nCC;
end
% 全部QDCT低频系数索引,一行代表一个低频系数的索引
indAllCoefLF = zeros(nBlock*nLF, 2);
offSet = (indBlock - ones(size(indBlock))) .* 8;
k=1;
for i=1:nBlock
    indAllCoefLF(k:k+nLF-1, :) = indSubBlockLF + ones(nLF,1) * offSet(i,:);
    k=k+nLF;
end

% 全部低频系数的个数
nAllLF = size(indAllCoefLF,1);
% 转换为线性索引
buf = zeros(nAllLF,1);
for i=1:nAllLF
    buf(i) = sub2ind([rQDCT,cQDCT], indAllCoefLF(i,1), indAllCoefLF(i,2));
end
indAllCoefLF = buf;
% 全部的低频系数
coeffLF = QDCT(indAllCoefLF);

% -------------------------------------------------------------------------
% 修改低频系数
% 高频系数索引
indCoefHF = setxor(indAllCoefLF, linspace(1,rQDCT*cQDCT,rQDCT*cQDCT));
% 可修改系数索引 indChangeable
mask = QDCT;
mask(indCoefHF) = inf;
% 可修改系数的索引
indZero = find(mask==0);
nZeroLF = length(indZero);
indChangeable = find(abs(mask) == 1);

% 可修改系数的个数
nChangeable = length(indChangeable);
% 被修改的系数个数
nChange = ceil(changeRate * nChangeable);
% 置乱
indChangeable = indChangeable(randperm(nChangeable));
indBeChanged = indChangeable(1:nChange);
% 修改绝对值为1的系数
QDCT(indBeChanged) = QDCT(indBeChanged) + sign(QDCT(indBeChanged));

% 修改零系数
if(~isnumeric(changeRateZero) || ~exist('changeRateZero','var'))
    changeRateZero=0;
else
    indZero = indZero(randperm(nZeroLF));
    nChangeZero = round(changeRateZero * nZeroLF);
    indBeChangedZero = indZero(1:nChangeZero);
    half = round(length(indBeChangedZero)/2);
    indN = indBeChangedZero(1:half);
    indP = indBeChangedZero(half+1:end);
    QDCT(indP) = 1;
    QDCT(indN) = -1;
end

% 逆向嵌入F5
QDCT = ReverseF5(QDCT, 0.0);
if(exist('debug', 'var'))
    fprintf('可修改系数1的个数%d\n',nChangeable);
    fprintf('低频零系数个数%d\n', nZeroLF)
    fprintf('被修改系数1的个数%d\n', nChange);
    fprintf('被修改系数0的个数%d\n', nChangeZero);
end
end

% 逆向嵌入F5
function QDCT = ReverseF5(QDCT, changeRate)
% 
% [QDCT, QTable, Cb, Cr] = getQDCT(cover, quality);
%% ----------------------------------------------------------------
% 非零AC系数
nzAC = nnz(QDCT)-nnz(QDCT(1:8:end,1:8:end));
% 被修改的系数的个数
nChanges = ceil(changeRate*nzAC);
if(nChanges<1)
    return;
end
% nChanges = ceil(nzAC / 2^k);
changeable = (QDCT~=0);
changeable(1:8:end,1:8:end) = false;
% 可修改的系数的位置
changeable = find(changeable);
%rng('state',seed);
% 置乱
changeable = changeable(randperm(nzAC)); % create a pseudorandom walk over nonzero AC coefficients
to_be_changed = changeable(1:nChanges);   % coefficients to be changed

% 绝对值加一    
QDCT(to_be_changed) = QDCT(to_be_changed) + sign(QDCT(to_be_changed));

end

function k = getEfficiency(payload)
% 计算嵌入效率k
switch(payload)
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
        fprintf('Unexpected payload:%f.2\n', payload);
        pause;
end
end