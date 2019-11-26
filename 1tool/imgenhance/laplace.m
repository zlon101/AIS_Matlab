function laplace()
% 拉普拉变换，用于图像锐化
%% 
% close all; clc;
c = 1;
Root = 'D:\MATLAB_Software\myInstall\bin\images\tmp\';
imgData=imread([Root, '7.bmp']);
outImg = [Root, '7_laplaced.bmp'];
srcImg=double(imgData);
[R, C]=size(srcImg);

%% laplace    
% % Modle=[0, 1, 0; 1, -4, 1; 0, 1, 0] .* -1;
Modle=[-1, -1, -1; -1, 8, -1; -1, -1, -1];
laplaceResult=zeros(R, C);
dstImg=double(zeros(R,C));
for i=2:R-1
    for j=2:C-1
        localBlock = [srcImg(i-1,j-1),srcImg(i-1,j),srcImg(i-1,j+1);
                      srcImg(i,j-1),srcImg(i,j),srcImg(i,j+1);
                      srcImg(i+1,j-1),srcImg(i+1,j),srcImg(i+1,j+1);];
		laplaceResult(i,j) = sum( sum(localBlock .* Modle) );        
        dstImg(i,j) = round( srcImg(i,j) + c*laplaceResult(i,j) );
    end
end
maxLaplaceGray=max(max(dstImg));
minLaplaceGray=min(min(dstImg));
dstImg=( (dstImg-minLaplaceGray) .* 255) ./ (maxLaplaceGray-minLaplaceGray);
dstImg=uint8(dstImg);

figure('name', 'dstImg');    imshow(dstImg);
figure('name', '高频部分');  imshow(laplaceResult);
% imwrite(dstImg, outImg, 'bmp');
end