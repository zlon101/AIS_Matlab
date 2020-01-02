function F = LoadFeature(FeatRoot)
% 加载从C++程序中提取出来的单幅图像的SRM特征
% 输出：符合检测器的特征数据
%% 
% FeatDir  = dir([FeatRoot '*.*']);         % 遍历所有**格式文件
% FeatDir(1)=[];   FeatDir(1)=[];
load('D:\Matlab2017b\bin\feature_extract\SRMFeatOrder.mat');% SRMFeatOrder
f = [];
% SRMFeatOrder = 1;
num = size(SRMFeatOrder, 1);
old='';
for i = 1:num                % 遍历结构体就可以一一处理图片了    
  path = [FeatRoot, SRMFeatOrder{i}, '.fea'];
  featmp = load(path);
  featmp(:, end) = [];
  f = [f, featmp];

  msg=sprintf('- count: %3d/%d', i, num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
% 图像的名字
path = [FeatRoot, SRMFeatOrder{1}, '.fea'];
featmps = load(path);  featmps(:, 1:end-1) = [];
names = cell(length(featmps), 1);
for i=1:size(names,1)
  names{i} = [num2str(featmps(i)), '.pgm'];
end
F.names = names;  F.F = f;
end