clc;close all;
FeatureDatas = 'D:\Program Files\Matlab2017b\bin\FeatureDatas\';
FcPath = [FeatureDatas,'\covers\','FC_BOSS_SRM.mat'];
FsPath = [FeatureDatas,'\stegos\','FS_HUGO_04_SRM.mat'];
MFc = matfile(FcPath);
MFs = matfile(FsPath);

% learner = TrainArgs(40);
FLD_Test(MFc, MFs);

% fprintf('S的投影中心：%.3f\n',mean(PS,1) * learner.w);
% fprintf('均值向量的距离:%.3f\n',norm(mean(PS)-mean(PC)));
% figure;
% learner = TrainArgs{minInd};
% PC = C_BOSS_SRM.F(:,learner.subspace) * learner.w;
% PS = S_HUGO_04_SRM.F(:,learner.subspace) * learner.w;
% border = learner.b;
% vmin = min(min(PC),min(PS));    vmax = max(max(PC),max(PS));
% axis([0,length(PC), vmin, vmax]);
% plot(PC, 1:size(PC,1), '.k');hold on;
% plot(PS, 1:size(PS,1), '.r');hold on;
% hp = plot([border, border],[0,size(PS,1)], '-r');hold on;
% legend('c','s', ['b: ', num2str(border)]);