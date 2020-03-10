function maxSRMEXE(imgRoot,outPath,startInd,endInd)
% 计算图像特征
% fetuStruct    全部维数特征
% addpath(genpath('./featureFile'));    % 添加文件路径
%%
imgDirs = dir([imgRoot, '*.pgm']);
num= length(imgDirs);
if(exist('startInd','var') && str2double(startInd)>0)
  startInd = single(str2double(startInd));
else
  startInd = 1;
end
if(exist('endInd','var'))
  endInd=single(str2double(endInd));
else 
  endInd=single(num);
end
names = cell(num,1);
for i=1:num
  names{i}=imgDirs(i).name;
end
names= sort(names); names=names(startInd:endInd);
clear imgDirs;
fprintf('# ind: %d - %d\n', startInd,endInd);

t0= datetime('now'); old='';
load('D:\Matlab2017b\bin\feature_extract\featureFile\maxSRM\maxSRMPropertys.mat',...
  'propertys');
% 一列为一个样本
F= zeros(34671, endInd-startInd+1,'single'); 
count= 1;
for i= startInd: endInd
  imgPath= [imgRoot,names{i}];
  I= imread(imgPath);
  F(:,count)= mergeProper(maxSRM(I,ones(size(I))), propertys);
  count=count+1;
  % 打印
  msg=sprintf('- count: %3d/%d',i,num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
fetuStruct.names= names;
fetuStruct.F= F';
save(outPath, 'fetuStruct');

fprintf('\nnumbel of img: %d\n', count-1);
fprintf('\n耗时: '); disp(datetime('now')-t0);
end

function F=mergeProper(Stru,propertys)
% propertys=fieldnames(Stru);
len=length(propertys);
F= zeros(34671,1,'single');
ind= 1;
for i=1:len
  T= getfield(Stru,propertys{i});
  F(ind:ind+length(T)-1)= T;
  ind= ind+length(T);
end
end