function laplace_zhihu()
close all; clc;
c = -1;
Root = 'D:\MATLAB_Software\myInstall\bin\images\tmp\';
imgData=imread([Root, '7.bmp']);
outImg = [Root, '7_laplaced.bmp'];
srcImg=double(imgData);
gray=srcImg;
[grayRow, grayColumn]=size(gray);

%% laplace
% laplacePlate=[0, 1, 0; 1, -4, 1; 0, 1, 0];
laplacePlate=[1, 1, 1; 1, -8, 1; 1, 1, 1];
laplaceResult=zeros(grayRow, grayColumn);
laplaceGray=zeros(grayRow, grayColumn);
laplaceGray=double(laplaceGray);
for i=1:grayRow
    for j=1:grayColumn
        for k=-1:1
            for n=-1:1
                if (i+k>=1) && (i+k<=grayRow) && (j+n>=1) && (j+n<=grayColumn)
                    grayValue=gray(i+k, j+n);
                else
                    grayValue=0;
                end
                laplaceResult(i, j)=laplaceResult(i, j)+laplacePlate(k+2, n+2)*grayValue;
            end
        end
        laplaceGray(i, j)=round(gray(i, j)+c*laplaceResult(i, j));
    end
end

maxLaplaceGray=max(max(laplaceGray));
minLaplaceGray=min(min(laplaceGray));
laplaceGray=((laplaceGray-minLaplaceGray) .* 255) ./(maxLaplaceGray-minLaplaceGray);
laplaceGray=uint8(laplaceGray);
fileName='../output/laplace';
fileSuf='.jpg';
gammaStr=num2str(c);
file=[fileName, gammaStr, fileSuf];

figure('name', 'gray');imshow(laplaceGray);
figure('name', 'laplace');imshow(laplaceResult);

% imwrite(laplaceGray, file);  
% imwrite(laplaceResult, '../output/laplaceResult.jpg');  
end