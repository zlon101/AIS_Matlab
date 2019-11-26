function [FcPcaed, FsPcaed, k] = PcaFunction(Fc,Fs,setDimension)
% PCA降维
% Fc、Fs     载体载密特征
% k         降维后的维度
% [coeff,score,latent,tsquared,explained,mu] = pca(X, Name, Value)
%    X: n*m 矩阵, 表示n个样本观测值, 每个样本为m维矢量, 要求n>=m
%    coeff: m*m的转换矩阵
%    score: 减维后的n*m矩阵, 一行对应一个样本, 一列对应一个主成分
%%
T = 0.98;
%% 标准化
[Fc, mu, sigma] = zscore(Fc);
Fs = bsxfun(@minus, Fs, mu);
Fs = bsxfun(@rdivide, Fs, sigma);

%% pca减维
[Tc,FcPcaed,latent_c,~,~,mu] = pca(Fc);
rank_c = cumsum(latent_c)./sum(latent_c);
FsTmp= bsxfun(@minus, Fs, mu);
FsPcaed = FsTmp * Tc;

% 两种减维方式, 那种更好
% border = size(Fc,1);
% [Tc,pcaed,latent_c,~,~,mu] = pca([Fc;Fs]);
% rank_c = cumsum(latent_c)./sum(latent_c);
% FcPcaed = pcaed(1:border,:);
% FsPcaed = pcaed(border+1:end,:);

%% 带权重的减维
%{
W = 1./var(Fc);
[Tc,FcPcaed,latent_c,~,~,mu] = pca(Fc,'VariableWeights',W);
rank_c = cumsum(latent_c)./sum(latent_c);
coefforth = diag(sqrt(W)) * Tc;
FsPcaed = zscore(Fs) * coefforth;
%}

for k=1:length(rank_c)
   if rank_c(k)>=T
       break;
   end
end
if(nargin<3)
    FcPcaed=FcPcaed(:,1:k);
    FsPcaed=FsPcaed(:,1:k);
else
    FcPcaed=FcPcaed(:,1:setDimension);
    FsPcaed=FsPcaed(:,1:setDimension);
end
end

%% 画图
function draw(Fc, Fs, dim, nSample)
% 画图:2维
if(dim==2)
    figure;
    scatter(Fc(1:nSample, 1), Fc(1:nSample, 2), 20, 'x', 'filled',...
        'MarkerEdgeColor','k', 'MarkerFaceColor','blue');
    hold on;
    scatter(Fs(1:nSample, 1), Fs(1:nSample, 2), 20, 'o', 'filled',...
        'MarkerEdgeColor','k', 'MarkerFaceColor','red');
    xlabel('X');ylabel('Y');zlabel('Z');
% 画图：3维
else
    hs1 = figure;
%     hs1 = subplot(2,1,1); hs2 = subplot(2,1,2);
    scatter3(Fc(1:nSample, 1), Fc(1:nSample, 2), Fc(1:nSample, 3), 20, 'x', 'filled',...
        'MarkerEdgeColor','k', 'MarkerFaceColor','blue');
    hold on;
    scatter3(Fs(1:nSample, 1), Fs(1:nSample,2), Fs(1:nSample,3), 20, 'o', 'filled',...
        'MarkerEdgeColor','k', 'MarkerFaceColor','red');
    xlabel('X');ylabel('Y');zlabel('Z');
    view(-0.7,90);       % 改变坐标轴的角度
end
%{
figure;
ind = 1;
for dim=1:reserve    
    subplot(1,2,ind);
    plot(FcPcaed(1:200, dim), '.r');hold on;
    plot(FsPcaed(1:200, dim), 'ob');
    title(['Fc & Fs in dimension ',num2str(dim)]);
    ind = ind + 1;
    if(mod(dim,2)==0)
        figure;
        ind = 1;
    end
end
%}
end