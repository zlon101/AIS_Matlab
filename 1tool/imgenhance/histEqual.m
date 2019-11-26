close all;
Root = 'D:\MATLAB_Software\myInstall\bin\images\tmp\';
coverPath = [Root, '7.pgm'];
% stegoPath = [Root, '7-HUGO05.pgm'];
% equaledCover = [Root, '7-Equaled.pgm'];
% equaledStego = [Root, '7_Equaled-HUGO05.pgm'];

% cover = imread(coverPath);
% stego = imread(stegoPath);

% coverEqual = histeq((cover));     % 直方图均衡化
% imwrite(coverEqual, [Root, 'coverEqualed.pgm'], 'pgm');

%% 提取特征

F = SRMexample({coverPath});

%%
%{
function histEqual()
f=imread('zftjhh1.jpg');
[m,n,d]=size(f);%灰度图1维，彩色图3维
if d==1
    f1=f;%复制后新的图片f1，作为改变后的图片
elseif d==3
    f=rgb2gray(f);
    f1=f;
end
figure 
imhist(f)
[count,x]=imhist(f);%count表示每个灰度级别有多少个像素，x表示有多少个灰度级别

PDF=count/(m*n);%PDF表示每个灰度级别出现的概率，一共有256行
CDF=cumsum(PDF);%CDF表示逐行相加的概率，也就是累加概率

for i=1:256
    xiangsuxushu=find(f==i);%原本灰度级别为i的像素在第几位
    changdu=length(xiangsuxushu);
    for j=1:changdu
        f1(xiangsuxushu(j))=round(CDF(i)*256-1);%每一个原本灰度级别为i的像素，
                                              %灰度级别改为累加出现概率*256
                                              %再取整
    end
end

figure
imhist(f1)
figure
imshow(f1)
end
%}