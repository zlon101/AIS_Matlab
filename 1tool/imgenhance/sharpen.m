function [imgData,HF] = sharpen(imgData, Amplitude)
% 图像锐化, 丰富纹理区域
% Amplitude: 幅度
%%
if(ischar(imgData))
   imgData = imread(imgData);
end
imgData = single(imgData);
[t, ~] = wiener2(imgData, [3, 3]);
HF = imgData - t;
% 限定修改位置只在纹理区域
% HF( abs(HF)<=T ) = 0;  HF(HF<0)=-1; HF(HF>0)=1; D = floor(HF * Amplitude);
T = 10;
HF = HF * Amplitude;  
% HF(HF>T) = T; HF(HF<-1*T) = -1*T;
imgData = (imgData + HF);
imgData(imgData<0) = 0;  imgData(imgData>255) = 255;