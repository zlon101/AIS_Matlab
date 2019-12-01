function [dstImg,HF] =  laplace(imgData, Am)
% 拉普拉斯变化, 用于图像锐化
% imgData: single
%% 
% % Modle=[0, 1, 0; 1, -4, 1; 0, 1, 0] .* -1;
Modle=[-1, -1, -1; -1, 8, -1; -1, -1, -1];
HF = filter2(Modle,imgData,'same');
% HF( abs(HF)<4 ) = 0;
HF = single(0.05*Am .* HF);
% T = 10; HF(HF>T) = T; HF(HF<-1*T) = -1*T;
dstImg = (imgData + HF);
dstImg(dstImg<0) = 0;  dstImg(dstImg>255) = 255;
end

% HF=zeros(R, C);
% for i=2:R-1
%     for j=2:C-1
%         localBlock = [imgData(i-1,j-1),imgData(i-1,j),imgData(i-1,j+1);
%                       imgData(i,j-1),imgData(i,j),imgData(i,j+1);
%                       imgData(i+1,j-1),imgData(i+1,j),imgData(i+1,j+1);];
% 		HF(i,j) = sum( sum(localBlock .* Modle) );
%     end
% end