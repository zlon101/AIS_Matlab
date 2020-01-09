function F = LoadFeature(FeatRoot)
% 加载从C++程序中提取出来的单幅图像的SRM特征
% 输出：符合检测器的特征数据
%% 
dimF=34671;  format='.pgm';

S=load('D:\Matlab2017b\bin\feature_extract\SRMFeatOrder.mat');% SRMFeatOrder
SRMFeatOrder=S.SRMFeatOrder;
ftmp = load( [FeatRoot,SRMFeatOrder{1},'.fea'] );
FT = zeros(size(ftmp,1), dimF, 'single');
num = size(SRMFeatOrder, 1);
old=''; Ind=1;
for i = 1:num
  path = [FeatRoot, SRMFeatOrder{i}, '.fea'];
  ftmp = load(path);
  ftmp(:, end) = [];
  FT(:, Ind:Ind+size(ftmp,2)-1) = ftmp;
  Ind = Ind+size(ftmp,2);
  % FT = [FT, ftmp];

  msg=sprintf('- count: %3d/%d', i, num);
  fprintf([repmat('\b',1,length(old)),msg]);  old=msg;
end
% 图像的名字
path = [FeatRoot, SRMFeatOrder{1}, '.fea'];
featmps = load(path);  featmps(:, 1:end-1) = [];
names = cell(length(featmps), 1);
for i=1:size(names,1)
  names{i} = [num2str(featmps(i)), format];
end
F.names = names;  F.F = single(FT);
end