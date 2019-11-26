% 各类其他程序, 主要用于验证各类计算是否正确
%% -------------------------------------------------------------------------

% S_F5_CCPEV_04_75    S_F5_CCPEV_02_75  
% S_nsF5_CCPEV_04_75  S_nsF5_CCPEV_02_75
% S_JUNIWD_CCPEV_04_75
C = C_CCPEV.F;  
S4 = S_nsF5_CCPEV_02_75.F;
S2 = S_JUNIWD_CCPEV_04_75.F;


%% PCA
load cities
w = 1./var(C);
[coeff,score,latent,tsquared,explained,mu] = ...
        pca(C);
% ,'VariableWeights',w
% 验证pca
x = bsxfun(@minus, C, mu);
score2 = x * coeff;
diff2 = round(score2 - score,10);
% 有权重VariableWeights时
coefforth = diag(sqrt(w)) * coeff;
cscores = zscore(C) * coefforth;
diff = round(cscores - score,10);


%% cover-stego的样本均值向量的欧式距离
mC = mean(C);
mS4 = mean(S4);     mS2 = mean(S2);
md4 = norm(mS4-mC);
md2 = norm(mS2-mC);

%% 计算不同算法对应的每个cover-stego样本对的特征的欧式距离
%   验证特征的欧式距离小, Pe是否就小
t4 = S4-C;   t2 = S2-C;
d4 = zeros(size(t4,1),1);
d2 = zeros(size(t2,1),1);
for i=1:size(t4,1)
    d4(i) = norm(t4(i,:));
end
for i=1:size(t2,1)
    d2(i) = norm(t2(i,:));
end
validate = d4 > d2;
v = nnz(validate==0);
fprintf('%d\n', v);



%% 对目录进行特征提取、分类
% imgDir  = 'E:\astego\embedAlg\ste_nsF5_0.4\';
% name = '321.jpg';
% quality = getQuality([imgDir, name]);
% Fc = ccpev548([coverPath, name], quality);
% S_nsF5_04_75 = getFeature(imgDir, quality, -3);
% embedInRoot(coverPath, payLoad, stegoPath,1);
% tutorial(C_CCPEV,S_JUNIWD_04_75);
% 用ccpev对F5、nsF5和J-UNIWD进行分析, 测试适应度是否与误检率一致
% 适应度：自定义L、欧式距离、减维后的欧式距离
% analyze_feature(C_CCPEV.F, S_F5_CCPEV.F, 'F5');
[fitsJuniwd] = calcu_fitness(C_CCPEV.F, S_JUNIWD_04_75.F);


%% 频域纹理, QDCT系数复杂部分
Q = getQuality(coverPath);
[cJobj, Q, ~] = getQDCT(coverPath);  cQdct = cJobj.coef_arrays{1};
nBlocks = size(cQdct,1) * size(cQdct,2) / 64;
% 计算每个8x8块的方差
% varBlock = @(block_struct) var(reshape(block_struct.data, [],1));
nnzOfBlockHandle = @(truct_data) nnz(truct_data.data)-1;
nnzOfBlock = (blockproc(cQdct, [8,8], nnzOfBlockHandle));    

% ind = sub2ind(nnzOfBlock,size(nnzOfBlock,1),size(nnzOfBlock,2));
[sordNnzOfBlock, index] = sort(nnzOfBlock(:), 'descend');
[r,c] = ind2sub(size(nnzOfBlock),index(4085));
r = r-1;  c = c-1;
% 被选中的QDCT块
blockQdct = cQdct(r*8+1:r*8+8, c*8+1:c*8+8); 

%% 统计块中非零Ac系数的个数
[count_c,center_c]=hist(nnzOfBlock, max(nnzOfBlock(:))-min(nnzOfBlock(:))+1 );
fig=figure;
binWidth = 1;
hb=bar(center_c-0.5*binWidth, count_c, binWidth, 'FaceColor', ([173,217,226])/255,...
    'EdgeColor','black', 'LineWidth',1.2, 'LineStyle', '-');
I = double(cover);
%% 二维自适应降噪滤波
J = round(wiener2(I, [3,3]));
diff = J-I; tmp = diff(1:8,1:8);
vm = mean(tmp(:));
figure; imshow(diff, []);
imwrite(cover(:,:,1), coverPath, 'jpg', 'Quality', 75);



%--------------------------------------------------------------------------
%  各类函数
%--------------------------------------------------------------------------
quality = getQuality(coverPath);
[jobj, quality, ~] = getQDCT(coverPath);  coverQdct = jobj.coef_arrays{1};


%% 标准化  sigma = std(x);
[C, mu, sigma] = zscore(C);
S = bsxfun(@minus, S, mu);
S = bsxfun(@rdivide, S, sigma);


%% 保存Jobj
jobj.coef_arrays{1} = DCT;
jobj.optimize_coding = 1;
jpeg_write(jobj,STEGO);
imwrite(cover(:,:,1), coverPath, 'jpg', 'Quality', 75);

%% 隐写算法
nsf5_simulation(coverPath, stegoPath, payLoad, seed);
F5(coverPath, stegoPath, payLoad);
J_UNIWARD(coverPath, stegoPath, single(payLoad));

%% 随机流stream
stream = RandStream.getGlobalStream 
prevstream = RandStream.setGlobalStream(stream)

randstream = RandStream('mt19937ar','Seed',sum(100*clock));
seed_subspaces = rand(randstream);
seed_bootstrap = rand(randstream);

randstream.subspaces = RandStream('mt19937ar','Seed',seed_subspaces);
randstream.bootstrap = RandStream('mt19937ar','Seed',seed_bootstrap);

%% 计算运行时间
t0 = datetime('now');
t1 = datetime('now');  fprintf('嵌入耗时: '); disp(t1-t0);

%%
indSubBlockLF = [8,8; 
             8,7; 7,8;
            8,6;7,7;6,8;
           8,5;7,6;6,7;5,8;
          8,4;7,5;6,6;5,7;4,8;
         8,3;7,4;6,5;5,6;4,7;3,8;
        8,2;7,3;6,4;5,5;4,6;3,7;2,8;
       8,1;7,2;6,3;5,4;4,5;3,6;2,7;1,8];
