function [imgData, D] = multiScaleSharpen(imgData, Am)
% ��߶�ͼ��ϸ����ǿ
%% 
imgData = single(imgData);
sigma1 = 1.0; sigma2 = 2.0; sigma3 = 4.0;
Radius = 2; % �˲�����С
H1 = fspecial('gaussian', [Radius,Radius], sigma1);
H2 = fspecial('gaussian', [Radius*2-1,Radius*2-1], sigma2);
H3 = fspecial('gaussian', [Radius*4-1,Radius*4-1], sigma3);
% B1,B2,B3 �ֱ�������ģ��֮���ͼ��
B1= imfilter(imgData, H1, 'replicate');
B2= imfilter(imgData, H2, 'replicate');
B3= imfilter(imgData, H3, 'replicate');
% D1,D2,D3����������ε�ϸ��
D1=imgData-B1; D2=B1-B2; D3=B2-B3;
%%
w1=0.5; w2=0.5; w3=0.25;
D = (1 - w1.*sign(D1)).*D1 + w2*D2 + w3*D3;
D = D .* Am;
T = 10;
D(D>T) = T; D(D<-1*T) = -1*T;
imgData = D + imgData;
imgData = (imgData);
imgData(imgData<0) = 0;  imgData(imgData>255) = 255;