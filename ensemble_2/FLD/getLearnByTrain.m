function [learner]=getLearnByTrain(FC, FS, learner)
% 输入: 载体载密图像特征
% 输出: 集成分类器的参数-投影向量，边界
% 验证两类特征的可区分性
%% 
if(exist('learner','var'))
  FC = FC(:, learner.subspace);
  FS = FS(:, learner.subspace);
end

%% 标准化  sigma = std(x);
% [FC, mu, sigma] = zscore(FC, 1, 1);
% FS = bsxfun(@minus, FS, mu);
% FS = bsxfun(@rdivide, FS, sigma);
% CMuSigma.mu = mu; CMuSigma.sigma = sigma;

%% FLD
if(~exist('learner','var'))
  learner = FLD_Ensemble(FC, FS);
end

%% 测试免疫处理的效果
% load('learner.mat');    load('CMuSigma.mat');
%{
vC = getValProjected(coverPath, learner, mu, sigma);
vS = getValProjected(stegoPath, learner, mu, sigma);
vS2 = getValProjected([imgRoot,'8-Juniwd.jpg'], learner, mu, sigma);
vT = getValProjected(immunedStegoPath, learner, mu, sigma);

pos = round(size(PC,1)*0.5);
plot(vC, pos, '*b');hold on;
plot(vS, pos, '*r');hold on;
plot(vS2, pos, 'dr');hold on;
plot(vT, pos, 'xr');hold on;
%}

%% 测试样本
%{
plot(PC2, nTrain+1:nSamp, '*b');hold on;
plot(PS2, nTrain+1:nSamp, 'or');hold on;
hp = plot([border, border],[0,nSamp], '-r');
legend('c','s','c','s',['b: ', num2str(learner.b)]);
title(algName,'Interpreter','none');
% legend(hp,{['b: ', num2str(learner.b)]}, 'FontSize',12);
%}
end

%% 计算单个图像的特征投影后的值
function vProjected = getValProjected(imgPath, learner, mu, sigma)
Q = getQuality(imgPath);
F = ccpev548(imgPath, Q);

C = bsxfun(@minus, F, mu);  C = bsxfun(@rdivide, C, sigma);
vProjected = C * learner.w;
end

%% 计算类中心与投影向量的角度
%{
cosangle = (mS * learner.w) / ( norm(mS)*norm(learner.w) );
angle = acos(cosangle) * 180 / pi;          % 角度
fprintf('diag: %f\n',angle);
%}