function ImageEnhance(inRoot,outRoot, Am)
%{
close all;clc;
bossDir='E:\astego\Images\BOSS_ALL\';
cRoot='E:\astego\StandExpers\covers\';
sRoot='E:\astego\StandExpers\HUGOlike_04\';
payload=single(0.3);

% names={'195.pgm';'103.pgm';'1004.pgm','1013.pgm'};
cDirs = dir([cRoot,'*.pgm']);
D0=zeros(length(cDirs),1);
psnr0=zeros(length(cDirs),1);
for i=1:length(cDirs)
  cover= single(imread([cRoot,cDirs(i).name]));
  stego= HUGO_like(uint8(cover), payload);stego=single(stego);
  D0(i)=calcuDist(cover,stego);
  psnr0(i)= cacul_psnr(cover, stego);
end
%}

%% ------------------------------------------------------------------
%{
cDirs = dir([cRoot,'*.pgm']);
names = cell(length(cDirs),1);
for i=1:length(cDirs)
  names{i} = cDirs(i).name;
end
D08 = zeros(length(names),1);
psnr08 = zeros(length(names),1);
for i=1:length(names)
  cover = single( imread([cRoot,names{i}]) );
  % 图像增强
  % [immuImg2,HF] = imgLaplace(coverImg, 1.5);
  [immuImg,HF] = sharpen(cover, 0.8);
  stego = HUGO_like(uint8(immuImg),payload); stego=single(stego);
  D08(i) =  calcuDist(immuImg,stego);
  psnr08(i)= cacul_psnr(cover, stego);
end

clear bossDir cDirs sDirs cRoot sRoot immuImg cover i payload stego ...
  names;
%}

%% 批量增强
% inRoot = 'E:\astego\Images\BOSS_ALL\';
% outRoot= 'F:\astego\Images\sharped\锐化_S31_Am0.8\';
Am = str2double(Am);

format = 'pgm';
dirs = dir([inRoot, '*.',format]);
names = cell(length(dirs),1);
old='';  total = length(dirs);
for i = 1:total
  names{i} = dirs(i).name;
  cpath=[inRoot, dirs(i).name];
  coverImg = single(imread(cpath));

  % sharpedData =  imgLaplace(coverImg, bestAbs{i,2});
  sharpedImg = sharpen(coverImg, Am);
  imwrite(uint8(sharpedImg), [outRoot,names{i}],format);

  msg=sprintf('- count: %3d/%d',i, total);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end