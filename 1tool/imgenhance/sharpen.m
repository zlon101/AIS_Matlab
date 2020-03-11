function [imgData,HF] = sharpen(imgData, Amplitude)
% 图像锐化, 丰富纹理区域
% Amplitude: 幅度
%%
if(ischar(imgData))
   imgData= single(imread(imgData));
end
[t, ~]= wiener2(imgData, [3,3]);
HF= imgData - t;
% 原T=8; T=3时, PSNR较好, 第一个算法时取消注释
T= 10; HF(HF>T)=T; HF(HF<-1*T)=-1*T;
HF= round(HF * Amplitude);
imgData= imgData + HF;
imgData(imgData<0)= 0;  imgData(imgData>255)= 255;