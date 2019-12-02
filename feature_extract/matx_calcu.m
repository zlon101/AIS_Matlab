function [T,Tais,constant]=matx_calcu(coeff,fc)
% 计算由cover-img-->steg-img的转换矩阵T，单维线性y=ax+b
% coeff:    一次项系数和常数项
% fc:       载体特征向量
% B:        常数项
% T：       n*n转换矩阵,n是特征维数
% T2:       n*n免疫矩阵
% a=coeff(:,1)   b=coeff(:,2);次数降序排列

ndw=size(coeff,1);      %特征维数
nsmp=size(fc,1);        %样本个数

% for i=1:ndw    
%     T(j,i)=a(i);        % 隐写算法矩阵    
% end
T=coeff(1:end-1,1:end);      %regress拟合得到
B=coeff(end,1:end);          % 1*ndw常数项

% 求免疫处理矩阵:AX=Y 则:X=A\Y
Y=[];
for i=1:nsmp
    Y=[Y;(fc(i,:)-B)];
end
X=fc\Y;
Tais=(T'\X')';
constant=B;             %常数项
end