function fetuStruct = getFeatures(imgRoot, nSamp)
% 计算图像特征
% imgPath       图像目录
% numSample:    设置样本个数
% fetuStruct    全部维数特征
% addpath(genpath('./featureFile'));    % 添加文件路径
%%
F=[];
count=0;  old='';
imgDirs = dir([imgRoot, '*.pgm']);  % 遍历所有**格式文件
% ImgLists(1)=[];ImgLists(1)=[];
if isempty(imgDirs)
  names = split(imgRoot, '\');
  names = names(end);
  F = structProce(SRM({imgRoot}),0);
else
  names = cell(length(imgDirs),1); % 图像名,不含全部路径
  if nSamp<=0
    nSamp=length(imgDirs);
  end
  for i = 1:length(imgDirs)
    imgPath=[imgRoot imgDirs(i).name];
    names{i}=imgDirs(i).name;
    F=[F;SRMProces( SRM({imgPath}), 0)];
    %F=[F;ccpev548(imgName,Q)];
    %F=[F;chen486(imgName)'];
    %F=[F;cchen972(imgName,80)'];
    %F=[F;SPAM686(imgName)'];
    %F=[F;mainDctrMex(imgName)];
    %F=[F;mainDctrMatlab(imgName)];
    %F=[F;CSR(imgName)];
    %F=[F;structProce(CFstar(imgName,80),1)];
    %F=[F;structProce(PSRM(imgName), 0)];
    
    % 打印
    count=count+1;
    msg=sprintf('- count: %3d/%d',count,nSamp);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
    if count>=nSamp
      break;
    end
  end
  fprintf('\nnumbel of img: %d\n', count);
end
[fetuStruct.names,ind]=sort(names);
fetuStruct.F=F(ind, :);
end