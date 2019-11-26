function FLD_Test(MFc, MFs, learner)
% 输入: 载体载密图像特征
% 输出: 集成分类器的参数-投影向量，边界
% 验证两类特征的可区分性
%% 
nTrain = 700;
C = MFc.FC_BOSS_SRM;
S = MFs.FS_HUGO_04_SRM;
if(exist('learner','var'))
    C = C(:, learner.subspace);
    S = S(:, learner.subspace);
end
% 两类的均值向量
mC = mean(C);   mS = mean(S);
D = norm(mC-mS);

%% 标准化  sigma = std(x);
[C, mu, sigma] = zscore(C, 1, 1);
S = bsxfun(@minus, S, mu);
S = bsxfun(@rdivide, S, sigma);
CMuSigma.mu = mu; CMuSigma.sigma = sigma;
% save('CMuSigma', 'CMuSigma');

%% FLD
nTrain = size(C, 1);
CTrain = C(1:nTrain,:);  STrain = S(1:nTrain, :);
%  S = C(nTrain+1:end,:);  STest = S(nTrain+1:end,:);
if(~exist('learner','var'))
    learner = FLD_Ensemble(CTrain, STrain);
end
%  标准化投影向量w
scale = 10/norm(learner.w);
% learner.b = learner.b .* scale;  % 保留2位小数
% learner.w = learner.w .* scale;  % 使w长度为10
fprintf('S的投影中心：%.3f\n',mean(S,1) * learner.w);
fprintf('均值向量的距离:%.3f\n',norm(mean(S)-mean(C)));
save('learner','learner');
%}

%% 可视化
%{
figure;
PC = C * learner.w;     PS = S * learner.w;     border = learner.b;
vmin = min(min(PC),min(PS));    vmax = max(max(PC),max(PS));
% axis([0,length(PC), vmin, vmax]);
plot(PC, 1:size(PC,1), '.k');hold on;
plot(PS, 1:size(PS,1), '.r');hold on;
hp = plot([border, border],[0,size(PS,1)], '-r');hold on;
legend('c','s', ['b: ', num2str(border)]);
axis([-65, -50, 0,length(PC)]);
%}

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