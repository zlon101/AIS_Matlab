function featuExtrct(inDir,varargin)
% F:\锐化_Am1.0_HUGO_03\ 1\  2\  E:\featData\**.mat
% 提取从C++程序中提取出来的特征数据
%%
dimF=34671;  nsamp=10000;
outPath = varargin{end};
if(nargin==2)
  Feat = LoadFeature(inDir);
else
  F=zeros(nsamp,dimF,'single'); names=cell(nsamp,1);
  Ind=1;
  for i=1:length(varargin)-1
    dir1=[inDir,varargin{i}];
    tmp = LoadFeature(dir1);
    F(Ind:Ind+size(tmp.F,1)-1, :) = tmp.F;
    names(Ind : Ind+size(tmp.F,1)-1) = tmp.names;
    Ind=Ind+size(tmp.F,1);
  end
  Feat.F=F; Feat.names=names;
  clear F names;
end
Feat.F=single(Feat.F);

% path1 = 'F:\astego\特征CPP\SUNWD_0.1\';
% F1 = LoadFeature(path1);
% Feat.names= F.names;  Feat.F= single(F.F);

[Feat.names, Ind]= sort(Feat.names);
Feat.F = Feat.F(Ind, :);
save(outPath, 'Feat');
fprintf('\n# end!');