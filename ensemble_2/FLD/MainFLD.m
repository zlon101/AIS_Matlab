
FeatureDatas = 'D:\Program Files\Matlab2017b\bin\FeatureDatas\';
FcPath = [FeatureDatas,'\covers\','C_BOSS_SRM.mat'];
FsPath = [FeatureDatas,'\stegos\SUNWD\','S_SUNWD_04_SRM.mat'];
% load(FcPath); FC = C_BOSS_SRM.F;
load(FsPath); FS = S_SUNWD_04_SRM.F;
clear C_BOSS_SRM S_SUNWD_04_SRM S_HILL_04_SRM S_HUGO_04_SRM...
  S_CZL2_04_SRM;

% [PE,bestLearner]=getBestLearner(learners, FC, FS);
PE = learnerTest(FC,FS,bestLearner);
% [learner]=getLearnByTrain(FC, FS);
% 缩放投影向量w
% scale = 10/norm(learner.w);
% learner.b = learner.b .* scale;  % 保留2位小数
% learner.w = learner.w .* scale;  % 使w长度为10
%}

DC = FC(:, bestLearner.subspace) * bestLearner.w;
DS = FS(:, bestLearner.subspace) * bestLearner.w;
D_SUNWD= DS-DC;
%% 可视化
%{
figure;
border = learner.b;
vmin = min([DC;DS]);  vmax = max([DC;DS]);
plot(DC, 1:size(DC,1), '.k');hold on;
plot(DS, 1:size(DS,1), '.r');hold on;
hp = plot([border, border],[0,size(DS,1)], '-b');hold on;
legend('c','s', ['b: ', num2str(border)]);
axis([vmin,vmax, 0,length(DC)]);
%}
clear vmin vmax hp FcPath FeatureDatas FsPath border ans;
%%