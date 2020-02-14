function z_test()
srcRoot = 'E:\astego\Images\standard_images\tif\';
dstRoot = 'E:\astego\Images\standard_images\bmp\';
numSample = -1;

imgFiles  = dir([srcRoot, '*.tif']);
Names=cell(length(imgFiles),1);
nImages = length(imgFiles);
if (~exist('numSample', 'var') || numSample<1)
  numSample=nImages;
end
old='';
for i = 1:nImages
  Names{i}=imgFiles(i).name;
  src = imread([srcRoot,Names{i}]);
  name = split(Names{i},'.');  name = name{1};
  imwrite(src(:,:,1), [dstRoot,name,'.bmp'], 'bmp');
  % 打印
  msg=sprintf('- count: %3d/%d',i,numSample);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;    
  if i>=numSample
      break;
  end
end
end

%% 频域隐写, nsF5
function ztest()
close all;
clear all;
seed = 99;rng('default');rng(seed);tic;   
imgRoot = 'E:\astego\embedAlg\img_jpg_500\';
imgName = '41';
payLoad = single(0.40);                          % bpac:非零Ac系数
coverPath = [imgRoot, imgName, '.jpg'];
payLoads = single([0.1,0.2,0.3,0.4]);

quality = getQuality(coverPath);
Fc = ccpev548(coverPath, quality);
%% 隐写
algNames = {'F5', 'nsF5', 'Juniwd'};
stegoPath = [imgRoot, 'z-',algNames{3},'-',imgName,'.jpg'];
algHands = {@F5, @nsf5_simulation, @J_UNIWARD};
getDiffOfFeature(coverPath,stegoPath,algHands{3}, payLoad)

Fs = cell(length(algHands), 1);
fits = zeros(length(algHands), 1);
% 不同算法
for i=1:length(algHands)    
    stegoPath = [imgRoot, 'z-',algNames{i},'-',imgName,'.jpg'];    
    algHands{i}(coverPath, stegoPath, payLoad);
    % 提取特征
    Fs{i} = ccpev548(stegoPath, quality);
    fits(i) = calcu_fitness(Fc, Fs{i});
    fprintf('%8s:  %.3f\n',algNames{i}, fits(i));
end
%}
% 不同嵌入率
%{
Fs = cell(length(payLoads), 1);
fits = zeros(length(payLoads), 1);
i = 1;
for j=1:length(payLoads)
    break;
    stegoPath = [imgRoot, 'z-',algNames{i},'-',imgName,'.jpg'];    
    algHands{i}(coverPath, stegoPath, payLoads(j));
    % 提取特征
    Fs{j} = ccpev548(stegoPath, quality);
    fits(j) = calcu_fitness(Fc, Fs{j});
    fprintf('%.2f:    %.3f\n',payLoads(j), fits(j));
end
%}
%% 提取特征

%% 特征分析
%{
nameAlg = 'F5'; nameFea = 'CCPEV';
Fc = ccpev548(coverPath, quality);
Fs = ccpev548(stegoPath, quality);
% Fts = ccpev548(tredStegoPath, quality);
analyze_feature(Fc, Fs,  '处理前');
% analyze_feature(Fc, Fts, '处理后');


Fc(166:168)=0;  Fs(166:168)=0;  Fts(166:168)=0;
Fc(166+274:168+274)=0;  Fs(166+274:168+274)=0;  Fts(166+274:168+274)=0;
% analyze_histogram(Fc(1:274), Fs(1:274), '处理前特征相似性');
% analyze_histogram(Fc(1:274), Fts(1:274),'处理后特征相似性');
analyze_histogram(Fc, Fs, '处理前特征相似性');
analyze_histogram(Fc, Fts,'处理后特征相似性');

% T = toc;    fprintf('耗时:%f\n', T);
%}

%% 直方图统计
%{
jobj = getQDCT(coverPath, quality);coverQDCT = jobj.coef_arrays{1};
jobj = getQDCT(stegoPath, quality);stegoQDCT = jobj.coef_arrays{1};
% jobj = getQDCT(treatedStegoPath, quality);treatedStegoQDCT = jobj.coef_arrays{1};
statistic(coverQDCT, stegoQDCT, -8:8);
% statistic(coverQDCT, treatedStegoQDCT, -8:8);

% 计算PSNR

ps_nr=cacul_psnr(coverName, stegoName)
ps_nr=cacul_psnr(coverName, treatedStegoName)
%}

%% 计算适应度
end

%% 对比载密与载体图像特征的相对变化量
function getDiffOfFeature(coverPath,stegoPath,algHandle,payload)
Q = getQuality(coverPath);
algHandle(coverPath, stegoPath, payload);
Fc = ccpev548(coverPath, Q);
Fs = ccpev548(coverPath, Q);
analyze_feature(Fc, Fs);
end