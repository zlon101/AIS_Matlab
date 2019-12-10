close all;clc;
%{
% coverRoot='E:\astego\Images\covers\锐化G1后载体_直接锐化版\';
% stegoRoot='E:\astego\Images\stegos\HUGO\锐化G1_HUGO_04_1万\';
names={'195.pgm';'103.pgm';'1004.pgm'};
bossRoot='E:\astego\Images\BOSS_ALL\';
% D1=zeros(length(names),1);
% psnr1=zeros(length(names),1);
% for i=1:length(names)
%   D1(i)= calcuDist(single(imread([coverRoot,names{i}])),...
%       single(imread([stegoRoot,names{i}])));
%   psnr1(i)= cacul_psnr(single(imread([bossRoot,names{i}])),...
%       single(imread([stegoRoot,names{i}])));
% end

%% ------------------------------------------------------------------
Root = 'E:\astego\Images\standard_test_images\bmp\';
dirs = dir([Root,'*.bmp']);
names = cell(length(dirs),1);
for i=1:length(dirs)
  names{i} = dirs(i).name;
end
D1 = zeros(length(names),1);
psnr1 = zeros(length(names),1);
for i=1:length(names)
coverImg = single( imread([Root,names{i}]) );
% 图像增强
% [immuImg2,HF] = imgLaplace(coverImg, 1.5);
[immuImg1,HF] = sharpen(coverImg, 1.1);
stegoImg1 = HUGO_like(uint8(immuImg1), single(0.4));stegoImg1=single(stegoImg1);
D1(i) =  calcuDist(immuImg1,stegoImg1);
psnr1(i) = cacul_psnr(coverImg, stegoImg1);
end
%}


%% 批量增强
inRoot = 'E:\astego\Images\BOSS_ALL\';
outRoot= 'E:\astego\Images\covers\锐化T2Cover\';
dirs = dir([inRoot, '*.pgm']); % imgDir(1)=[];imgDir(1)=[];
names = cell(length(dirs),1);
old='';
total = length(dirs);
for i = 1:total
  names{i} = dirs(i).name;
  path=[inRoot, dirs(i).name];
  coverImg = single(imread(path));

  % sharpedData =  imgLaplace(coverImg, bestAbs{i,2});
  sharpedImg = sharpen(coverImg, 1.1);
  imwrite(uint8(sharpedImg), [outRoot,names{i}], 'pgm');

  msg=sprintf('- count: %3d/%d',i, total);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
%}