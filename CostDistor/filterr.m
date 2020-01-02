%% 滤波器

H = fspecial('gaussian',3,0.5);  % gaussian
Y1 = imfilter(I,H);

%% SRM提取残差所用的高通滤波器
F = [-1, 2, -2, 2,-1;
      2,-6,  8,-6, 2;
     -2, 8,-12, 8,-2;
      2,-6,  8,-6, 2;
     -1, 2, -2, 2,-1];
   
%% HILL所用KB高通滤波器
KB= [-1, 2,-1;
      2,-4, 2;
     -1, 2,-1];
   
%% 维也纳低通滤波器
J = round(wiener2(I, [3,3]));

L= [1,2,1;
    2,4,2;
    1,2,1];
  
%% 高斯低通
% B = imgaussfilt(A)
montage(I)
%% 均值滤波