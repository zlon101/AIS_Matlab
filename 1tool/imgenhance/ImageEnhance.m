close all;clc;
% clear;clc;
Root = 'E:\astego\Images\test\';
srcImg = single(imread([Root, 'cover\195.pgm']));

D0= calcuDist(single(imread([Root, '锐化载体\195.pgm'])),...
    single(imread([Root, '锐化含密\195.pgm'])));
% psnr0= cacul_psnr(imread([Root, 'cover\195.pgm']),...
%     imread([Root, '锐化含密\195.pgm']));
%% 图像增强
[sharpImg1,HF1] = sharpen(srcImg, 1);
[sharpImg2,HF2] = imgLaplace(srcImg, 1);

% 隐写
stegoImg1 = HUGO_like(uint8(sharpImg1), single(0.4));stegoImg1=single(stegoImg1);
stegoImg2 = HUGO_like(uint8(sharpImg2), single(0.4));stegoImg2=single(stegoImg2);
% 性能计算
D1 =  calcuDist(sharpImg1,stegoImg1);
D2 =  calcuDist(sharpImg2,stegoImg2);
psnr1 = cacul_psnr(srcImg, stegoImg1);
psnr2 = cacul_psnr(srcImg, stegoImg2);
% figure('name','HF0'); imshow(round(HF0));
% figure('name','HF-1'); imshow(round(HF1));
% figure('name','HF2'); imshow(round(HF2));
%}

%% 批量增强
inRoot = 'E:\astego\Images\BOSS_ALL\';
outRoot = 'E:\astego\Images\covers\锐化G3载体_Am1\';
imgDir  = dir([inRoot '*.*']);    % 遍历所有**格式文件
imgDir(1)=[];imgDir(1)=[];
Names = cell(length(imgDir),1);    % 图像名,不含全部路径

old='';
total = length(imgDir);
for i = 1:total                % 遍历结构体就可以一一处理图片了    
    imgPath=[inRoot, imgDir(i).name];
    Names{i} = imgDir(i).name;
    srcImg = double(imread(imgPath));
    [sharpImg2, ~] = sharpen(srcImg, 1);
    imwrite(uint8(sharpImg2), [outRoot,Names{i}], 'pgm');
    
    msg=sprintf('- count: %3d/%d',i, total);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
end
%}