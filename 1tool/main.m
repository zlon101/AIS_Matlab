% HUGO隐写, 图像增强
%%
close all;
Root = 'E:\astego\Images\test\';
name = 'house.bmp';
payLoad = single(0.4);
srcPath = [Root, name];
srcStegoPath = [Root,'stegos\',name];
srcData = single(imread(srcPath));
% srcStegoData = HUGO_like(uint8(srcData), payLoad);
% imwrite(uint8(srcStegoData),srcStegoPath, 'pgm');
% 锐化
Am = 1.479;
[sharpedData, HF] = sharpen(srcData, Am);
% D= sharpedData-srcData; D(D==0)=NaN; figure;histogram(D);
% 隐写
sharpedStegoData = HUGO_like(uint8(sharpedData), payLoad);

sharpedPath = [Root,'sharpeds\',name];
sharpedStegoPath = [Root,'stegos\',name];
imwrite(uint8(sharpedData), sharpedPath, 'bmp');
imwrite(uint8(sharpedStegoData),sharpedStegoPath, 'bmp');

%% 提取特征
% fetuStruct = getFeatures(srcPath);  Fc = fetuStruct.F;
% fetuStruct = getFeatures(srcStegoPath);  Fs = fetuStruct.F;
% fetuStruct = getFeatures(sharpedPath);  Fc2 = fetuStruct.F;
% fetuStruct = getFeatures(sharpedStegoPath);  Fs2 = fetuStruct.F;
% fprintf('Am=%.2f  Fs2-Fc2=%.3f  Fs-Fc=%.3f\n', ...
%     Am,norm(Fs2- Fc2),norm(Fs- Fc));

%% PSNR
PSNR = cacul_psnr(sharpedStegoPath, srcPath); PSNR = round(PSNR,3);
%}