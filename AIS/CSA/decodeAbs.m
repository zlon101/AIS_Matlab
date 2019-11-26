function Abs = decodeAbs(genes,NumParas,Vmin,Vmax)
% genes:  编码抗体
% Val:    real value (precision: 6)
%%
[N, LS] = size(genes);
L = LS/NumParas;
Abs = [];
aux = 0:1:L-1;
aux = ones(N,1) * aux;    % 与Abs同大小的指数矩阵: [0,1,2;0,1,2]
for i=1:L:LS
    x = sum( (genes(:,i:i+L-1) .* 2.^aux), 2);    % 一行之和
    Vals = Vmin + x .* ( (Vmax-Vmin)/(2.^L-1) );
    Vals = round(single(Vals), 3);  % 列向量
    Abs = [Abs, Vals];  % 一行为一个样本
end
Abs = mat2cell(Abs, ones(size(Abs,1),1));  % N*NumParas cell
end