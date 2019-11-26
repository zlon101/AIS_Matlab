function [dstImg] = multiScaleSharpen(src, Radius)
% 多尺度图像细节增强
%% 
sigma1 = 1.0;
sigma2 = 2.0;
sigma3 = 4.0;
H1 = fspecial('gaussian', [Radius,Radius], sigma1);
H2 = fspecial('gaussian', [Radius*2-1,Radius*2-1], sigma2);
H3 = fspecial('gaussian', [Radius*4-1,Radius*4-1], sigma3);
% B1,B2,B3 分别是三次模糊之后的图像
B1= imfilter(src, H1, 'replicate');
B2= imfilter(src, H2, 'replicate');
B3= imfilter(src, H3, 'replicate');
% D1,D2,D3代表三个层次的细节
D1=src-B1;
D2=B1-B2;
D3=B2-B3;
%%
w1=0.5;
w2=0.5;
w3=0.25;
dstImg=(1 - w1.*sign(D1)).*D1 + w2*D2 + w3*D3 + src;
end