function featuExtrct(inDir,varargin)
% 提取目录中所有图像的特征
% F:\astego\特征CPP\锐化_Am1.0_HUGO_03\   2\  E:\featData\**.mat
%%
% cDir= 'E:\astego\Images\standard_images\covers\';
% sDir = 'E:\astego\Images\StandExpers\czl4\';
% name= '195.pgm';
% cPath=[sDir,name];
% F1= SRMProces(SRM(cpath), 0);
% S_CZL4_SRM_04 = getFeatures(sDir, -1);

%% 加载从C++程序中提取出来的特征数据
outPath = varargin{end};
if(nargin==2)
  Feat = LoadFeature(inDir); Feat.F=single(Feat.F);
else
  F=[]; names={};
  for i=1:length(varargin)-1
    dir1=[inDir,varargin{i}];
    tmp = LoadFeature(dir1);
    F=[F; tmp.F];
    names=[names; tmp.names];
  end
  Feat.F=F; Feat.names=names;
  clear F names;
end

% path1 = 'F:\astego\特征CPP\SUNWD_0.1\';
% F1 = LoadFeature(path1);
% Feat.names= F.names;  Feat.F= single(F.F);

[Feat.names, ind]= sort(Feat.names);
Feat.F = Feat.F(ind, :);
save(outPath, 'Feat');