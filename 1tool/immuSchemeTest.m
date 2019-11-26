function immuSchemeTest()
clc;
startTime = datetime('now');
%% 计算修改量的函数
calculNumOfModified = @(cQdct,sQdct) nnz((cQdct-sQdct));
%% 创建多个独立的随机流
[streamAlg,streamAntb,s3] = RandStream.create('mrg32k3a', 'NumStreams', 3);

imgRoot = 'D:\MATLAB_Software\myInstall\bin\images\tmp\';
imgName = '8';
payLoad = single(0.40);                          % bpac:非零Ac系数
coverPath = [imgRoot,imgName,'.jpg'];
%%
indAlg = 2;
algNames = {'F5', 'nsF5', 'Juniwd'};
algHands = {@F5, @nsf5_simulation, @J_UNIWARD};
stegoPath = [imgRoot,imgName,'-',algNames{indAlg},'.jpg'];
immunPath = [imgRoot,imgName,'-immuned','.jpg'];
immunedStegoPath = [imgRoot,imgName,'-immuned-',algNames{indAlg},'.jpg'];

[cJobj, Q, ~] = getQDCT(coverPath);  cQdct = cJobj.coef_arrays{1};
Fc = ccpev548(coverPath, Q);
algObj = struct('algOperat',algHands{indAlg},'Q',Q,'Fc',Fc,'payLoad',payLoad,...
    'cJobj',cJobj, 'immunedStegoPath', immunedStegoPath);
%{
algHands{indAlg}(coverPath, stegoPath, payLoad);
Fs = ccpev548(stegoPath, Q);
fprintf('未处理\n dist: %.3f\n', calcu_fitness(Fc, Fs));
% [~, dist] = analyze_feature(Fc, Fs,  '处理前');
sQdct = getQDCT(stegoPath); sQdct = sQdct.coef_arrays{1};
fprintf('cQdct-sQdct: %.d\n', calculNumOfModified(cQdct,sQdct));
fprintf('------------------------------\n');
%}
%%
[rImg, cImg] = size(cQdct);
nBlocks = size(cQdct,1) * size(cQdct,2) / 64;

%% 选中非零AC系数的个数超过threshold的8*8(块的二维索引)
%  计算每个8x8块中非零AC系数的个数
nnzOfBlockHandle = @(truct_data) nnz(truct_data.data)-1;
numnzOfBlock = (blockproc(cQdct, [8,8], nnzOfBlockHandle));    
[sordNnzOfBlocks, indexBlocks] = sort(numnzOfBlock(:), 'descend');
%  被选中的QDCT块**********************
% thresholdInd = round(0.1 * length(indexBlocks));
thresholdInd = 10;
indOfBlockSelected = indexBlocks(1:thresholdInd);
[r,c] = ind2sub(size(numnzOfBlock), indOfBlockSelected);
r = r-1;  c = c-1;
%  被选中的系数
selectedAllBlocksCoeff = zeros(size(cQdct));
indSelectedCoeff = [];
for i=1:length(r)
    t = zeros(size(cQdct));
    t( r(i)*8+1:r(i)*8+8, c(i)*8+1:c(i)*8+8) = 1;
    indSelectedCoeff = [indSelectedCoeff; find(t==1)];
    selectedAllBlocksCoeff( r(i)*8+1:r(i)*8+8, c(i)*8+1:c(i)*8+8) = 1;
end
%  被选中的块的全部系数的索引集合
%  indSelectedCoeff = find(selectedAllBlocksCoeff==1);
%  平滑块的系数的线性索引
indNoSelectedBlocksCoeff = find(selectedAllBlocksCoeff==0);

%% indOfExclude: 所有不符合条件的系数的线性索引
t1 = zeros(size(cQdct)); 
t1(1:8:end,1:8:end) = 1;
indOfDc = find(t1);                                 % DC系数的线性索引
indOfZeros = find(cQdct==0);                        % 零系数的线性索引
indOfExcludes = unique([indOfDc; indOfZeros; indNoSelectedBlocksCoeff]);
indSelectedCoeff = setdiff(indSelectedCoeff, indOfExcludes, 'stable');

%% 免疫处理方案

%  512*2矩阵，第一列表示列编码，第二列表示行编码，暂时只考虑行列相同的图像
%  零和DC系数的二维索引
%% 计算得到modifiedMap
memory = struct('antibs',[],'fits',[]);
pInit = 0.2;
nAntib = 15;                            % 抗体个数
nGener = 10;
bestAntibs = [];
bestFits = zeros(nGener,1);
nReplace = round(nAntib*0.3);
lenCode = length(indSelectedCoeff);
antibs = initCodes(nAntib, lenCode, pInit, streamAntb);
for i=1:nGener
    %% 计算一代种群的适应度       
    [fits,memory] = populaOperat(antibs, indSelectedCoeff, algObj, memory);
    [fits, ind] = sort(fits, 'ascend');     antibs = antibs(ind, :);
    %% 收集最佳基因片段
    priorGenes = gatherPriorGenes(antibs);
    
    %% 克隆 变异
    [pClones, pMus] = calcuPcloneMu(fits, pInit);        
    newAntibs = cloneMuOperat(antibs, priorGenes, pClones, pMus, streamAntb);   
    newAntibs = unique(newAntibs, 'rows');
    %% 第二阶段
    [newFits,memory] = populaOperat(newAntibs, indSelectedCoeff, algObj, memory);
    allFits = [fits; newFits];  allAntibs = [antibs; newAntibs];
    [allFits, indFits] = sort(allFits, 'ascend');
    allAntibs = allAntibs(indFits, :);
    priorGenes = gatherPriorGenes(allAntibs);
        
    bestAntibs = [bestAntibs; allAntibs(1,:)];
    bestFits(i) = min(allFits);
    %% 用新的抗体替换部分抗体
    allAntibs(nAntib+1:end, :) = [];    allFits(nAntib+1:end, :) = [];
    allAntibs(end-nReplace+1:end, :) = initCodes(nReplace, lenCode, pInit, streamAntb);

    [pClones, pMus] = calcuPcloneMu(allFits, pInit);    
    antibs = cloneMuOperat(allAntibs,priorGenes, pClones, pMus, streamAntb);
    antibs(end,:) = bestAntibs(end, :);
    antibs = unique(antibs, 'rows');
    fprintf('  %d:   %.3f\n', i, bestFits(i));
    fprintf('---------------------------------------------\n');
end
%}

%% 测试抗体
%{
load('bestAntibs.mat');
antib = bestAntibs(end, :);
fit = populaOperat(antib, indSelectedCoeff, algObj);
[D1, D2] = calcu_antibperfom(coverPath, stegoPath, immunedStegoPath);
analyze_feature(coverPath, [imgRoot,'8-Juniwd.jpg']);
%}

%% 计算运行时间
endTime = datetime('now');
duraTime = endTime-startTime;
% minuTime=minutes(duraTime);
fprintf('运行时长：');disp(duraTime);
%}
save('bestFits','bestFits');
save('bestAntibs','bestAntibs');
save('memory','memory');
end

%% 一次种群的计算过程
function [fits, memory] = populaOperat(antibs, indSelectedCoeff, algObj,memory)
cJobj = algObj.cJobj;   Fc = algObj.Fc;     Q = algObj.Q;
payLoad = algObj.payLoad;
immunedStegoPath = algObj.immunedStegoPath;
cQdct = cJobj.coef_arrays{1};
nAntib = size(antibs, 1);
fits = zeros(nAntib,1);
old_msg = '';
for j=1:nAntib
    %  除去已经计算过的抗体
    if(~isempty(memory.antibs))
        [~,ind] = ismember( antibs(j,:), memory.antibs, 'rows');
        if( ind>0 )
            fits(j) = memory.fits(ind);
            continue;
        end
    end
    immunQdct = cQdct;
    immunQdct(indSelectedCoeff) = immunQdct(indSelectedCoeff) + antibs(j,:)';
    immunJobj = cJobj;  immunJobj.coef_arrays{1} = immunQdct;

    %% 隐写
    algObj.algOperat(immunJobj, immunedStegoPath, payLoad);
    fits(j) = calcu_fit(Fc, immunedStegoPath, Q);
    %  记录抗体及其适应度
    memory.antibs = [memory.antibs; antibs(j,:)];
    memory.fits = [memory.fits; fits(j,:)];
    %  打印    
    msg = sprintf('- count: %d        %2.3f\n',j, fits(j));
    fprintf([repmat('\b',1,length(old_msg)), msg]);
    old_msg = msg;
end
fprintf('\n');
end

%% 一个抗体的计算过程
function fit = calcu_fit(Fc, immunedStegoPath, Q)
%% 图像质量
% PSNR = cacul_psnr(coverPath, immunedStegoPath);
% fprintf('PSNR: %.3f\n', PSNR);

%% 特征相似性
Fts = ccpev548(immunedStegoPath, Q);
fit = calcu_fitness(Fc, Fts);
end

%% 收集最佳基因片段
function [priorGenes,geneStatic] = gatherPriorGenes(antibs)
%% 方案1：比较适应度相邻两个抗体, 提取不同的部分

n = size(antibs,1);
geneStatic = zeros(n-1,size(antibs,2)) + nan;
for i=1:n-1
    ind = find(antibs(i,:) - antibs(i+1,:));
    geneStatic(i,ind) = antibs(i,ind);
end
[counts,centers] = hist(geneStatic, [-1:1]);
[~,ind] = max(counts);
priorGenes = centers(ind)';


%% 方案2：比较适应度相邻两个抗体, 提取相同的部分
%{
priorGenes = zeros(1,size(antibs,2)) + nan;
ind = find( (antibs(1,:)-antibs(2,:)) == 0 );
priorGenes(ind) = antibs(1,ind);
geneStatic =0;
%}
end

%% 每个抗体对应的克隆概率和变异概率
function [pClones, pMus] = calcuPcloneMu(fits, pInit)
bestFit = min(fits);
maxFit = max(fits);
range = max(fits) - min(fits);
nAntib = length(fits);
pClones = zeros(nAntib,1);
pMus = zeros(nAntib,1);
for i=1:nAntib
    pClones(i) = ((nAntib-i+1)/nAntib) * (maxFit - fits(i)) / range;
    pMus(i) = (fits(i) - bestFit) / range;
end
pClones = pClones * 0.5;
pClones(isnan(pClones)) = pInit;
pClones(isinf(pClones)) = pInit;

pMus(1) = pMus(2);
pMus(isnan(pMus)) = pInit;
pMus(isinf(pMus)) = pInit;
end

%% 克隆
function newAntibs = cloneMuOperat(antibs,priorGenes, pClones, pMus,randStrm)
num = size(antibs, 1);
newAntibs = [];
for i=1:num
    n = ceil(num*pClones(i));
    t = ones(n, 1) * antibs(i,:);
    t = muOperat(t, pMus(i), randStrm);
    % 注入最佳基因片段
    ind = find( isnan(priorGenes) );
    t(1:ceil(0.2*n), ind) = ones(ceil(0.2*n),1) * priorGenes(ind);
    newAntibs = [newAntibs; t];
end
newAntibs = unique(newAntibs, 'rows');
end

%% 变异
function antibs = muOperat(antibs, P, randStrm)
% pMu越大，变异的概率越大
num = size(antibs,1);
codes = single(rand(randStrm, num, size(antibs,2), 'single'));
codes(codes < P*0.5) = -1;
codes(codes > P) = 0;
codes(codes > 0) = 1;
antibs = mod(antibs + codes, 2);
end

%% 初始化抗体编码
function codes = initCodes(num, lenCode, P, randStrm)
% P越大，修改的概率越大
codes = [];
PS = [0.2, 0.4, 0.6, 0.8];
for i=1:4
    code = single(rand(randStrm, floor(num/4), lenCode, 'single'));
    code(code < PS(i)*0.5) = -1;
    code(code > PS(i) ) = 0;
    code(code > 0) = 1;
    codes = [codes; code];
end
n = num - size(codes,1);
if(n>0)
    code = single(rand(randStrm, n, lenCode, 'single'));
    code(code < P*0.5) = -1;
    code(code > P ) = 0;
    code(code > 0) = 1;
    codes = [codes; code];
end
end


%%  方案1：两列codes
%{
%% 有效编码索引
t2 = ones(size(cQdct))*111;
t2(indOfExcludes) = 0;
%  有效行编码
validRowsInd = find( sum(t2,2) > 0 );
%  有效列编码
validColsInd = find( sum(t2) > 0 );


%% 初始化编码
function [rCodes,cCodes] = initCodes(rows,cols,P, validRowsInd, validColsInd)
r = length(validRowsInd);
c = length(validColsInd);
validCodes = single(rand(r+c, 1,'single') > P);

rCodes = zeros(rows,1);
cCodes = zeros(cols,1);
rCodes(validRowsInd) = validCodes(1:r);    
cCodes(validColsInd) = validCodes(r+1:end);
end

%% 更新编码，变异
function [rCodes, cCodes] = updataCodes(rCodes, cCodes,P, validRowsInd,validColsInd)
rLen = length(validRowsInd);
cLen = length(validColsInd);
t = single(rand(rLen+cLen, 1,'single') < P);

rCodes(validRowsInd) = mod( rCodes(validRowsInd) + t(1:rLen), 2);
cCodes(validColsInd) = mod( cCodes(validColsInd) + t(rLen+1:end), 2);
end

%% 将编码转换为修改map
function modifiedMaps = codes2modfiedMap(rCodes, cCodes)
rLen = length(rCodes);
cLen = length(cCodes);
v = zeros(rLen, cLen);
for i=1:rLen
    for j=1:cLen
        v(i,j) = rCodes(i) * 2 + cCodes(j);
    end
end
modifiedMaps = rem(v-2, 2);
end
%}






