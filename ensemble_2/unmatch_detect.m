function [PE, PFA, PMD] = unmatch_detect(CTrn,STrn,CTst,STst)
% 非匹配检测, 载体库失配
%{
if(~isequal(CTrn,STrn))
  disp('错误');pause;
end
%}
%%
% load('trained_ensemble.mat');  load('results.mat');
% d_sub = results.optimal_d_sub;
settings= struct('verbose',2); % 'd_sub',d_sub,
num= 1;  PMD= zeros(num,1); PFA=zeros(num,1);
for seed = 1:num
  RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));
  random_permutation = randperm(size(CTrn,1));
  training_set = random_permutation(1:round(size(CTrn,1)/2));
  testing_set = random_permutation(round(size(CTrn,1)/2)+1:end);
  Tst_cover = CTst(testing_set,:);
  Tst_stego = STst(testing_set,:);
  [trained_ensemble,~]=ensemble_training(CTrn(training_set,:), STrn(training_set,:),settings);
  test_results_cover = ensemble_testing(Tst_cover,trained_ensemble);
  test_results_stego = ensemble_testing(Tst_stego,trained_ensemble);

  false_alarms = sum(test_results_cover.predictions~=-1);
  missed_detections = sum(test_results_stego.predictions~=+1);
  PFA(seed) = false_alarms / size(Tst_cover, 1);
  PMD(seed) = missed_detections / size(Tst_stego, 1);
  fprintf('PMD %i: %.4f\n',seed, PMD(seed));
end
PE = (PMD+PFA)/2;
disp('#10次测试的平均值:'); 
fprintf('PFA: %.4f+/-%.4f\n', mean(PFA),std(PFA));
fprintf('PMD: %.4f+/-%.4f\n', mean(PMD),std(PMD));
fprintf('PE:  %.4f+/-%.4f\n', mean(PE),std(PE));
