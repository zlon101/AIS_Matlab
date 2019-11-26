function [xaxis,yaxis]=DrawLineChart(D)
% 画整数统计折线图
% D:数据
D = D(:);
xaxis = min(D):max(D);
yaxis = zeros(length(xaxis),1);

j=1;
for v=xaxis(1):1:xaxis(end)
    yaxis(j) = length(find(D==v));
    j=j+1;
end

figure;
xlabel('Data');
ylabel('Frequency');
plot(xaxis,yaxis, '-kx')

% 频数最大的值和索引
% [maxFrequency, maxInd] = max(yaxis);
end