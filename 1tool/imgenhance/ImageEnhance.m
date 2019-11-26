% close all;
clear;clc;
Root = 'E:\astego\Images\test\';
imgData = double(imread([Root, '195.pgm']));
%% 图像增强
%方法1
% Radius=5;
% dstImg = multiScaleSharpen(src, Radius);
%方法2
% [dstImg, highFreq] = sharpen(imgData, 1);
% figure('name','sharpened');imshow(uint8(dstImg));
% figure('name','highFreq');imshow(highFreq, []);
% imwrite(uint8(dstImg), [Root, '195_Sharp.pgm'], 'pgm');
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
    imgData = double(imread(imgPath));
    [dstImg, ~] = sharpen(imgData, 1);
    imwrite(uint8(dstImg), [outRoot,Names{i}], 'pgm');
    
    msg=sprintf('- count: %3d/%d',i, total);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
end
%}