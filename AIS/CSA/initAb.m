function genes = initAb(num,len)
% 初始化抗体
% num: 抗体个数
% len: 抗体编码长度
%%
num = round(num);
genes = 2 .* rand(num,len) - 1;
genes(genes<0) = 0;
genes(genes>0) = 1;
% genes = hardlim(genes);  % 以零为阈值,映射为{0,1}
% Ab = hardlims(Ab);  % 以零为阈值,映射为{-1,1}
end