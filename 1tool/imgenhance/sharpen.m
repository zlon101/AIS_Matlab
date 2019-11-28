function [dstImg, HF] = sharpen(srcImg, Amplitude)
% 图像锐化，丰富纹理
% Amplitude: 幅度
%%
srcImg = single(imread(srcImg));
[t, ~] = wiener2(srcImg, [3, 3]);
HF = srcImg - t;
T = 2;
% 限定修改的位置只在纹理丰富的区域
% HF( abs(HF)<=T ) = 0;  HF(HF<0)=-1; HF(HF>0)=1; D = floor(HF * Amplitude);
T = 10;
D = HF * Amplitude;  
D(D>T) = T; D(D<-1*T) = -1*T;
dstImg =uint8(srcImg + D);
dstImg = single(dstImg);
% figure('name','HF');imshow(HF,[]);
% D(D==0)=NaN; figure('name', 'Diff'); histogram(D);
end