function [T,pcs] = reprod(genes, ClonedNum, fat)
% CloneNum-> 被克隆的数量
% fat	-> multiplying factor
% ind	-> best individuals
% T		-> temporary population
% pcs	-> final position of each clone
% cs：克隆cs倍
%%
T = [];
NumTotal = size(genes,1);
if ClonedNum == 1
    T = ones(NumTotal,1) * genes(1,:);
    pcs = NumTotal;
else
    cs = zeros(1,ClonedNum);
    pcs = zeros(1, ClonedNum);
    for i=1:ClonedNum
      % cs(i) = round(fat*N/i);      
      cs(i) = round(fat*NumTotal);
      pcs(i) = sum(cs);
      T = [T; ones(cs(i),1) * genes(i,:)];
    end
end
end