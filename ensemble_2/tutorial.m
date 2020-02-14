function [PE,PFA,PMD,trained_ensemble] = tutorial(cover, stego, numOfTrain)
% -------------------------------------------------------------------------
% names = intersect(cover.names,stego.names);    % 交集, 并去除重复项
% names = sort(names);
% 修改
names = cover.names;
% Prepare cover features C, imgCove与imgSteg要对应
cover_names = cover.names(ismember(cover.names,names));
[cover_names,ix] = sort(cover_names);
C = cover.F(ismember(cover.names,names),:);
C = C(ix,:);
% Prepare stego features S
stego_names = stego.names(ismember(stego.names,names));
[stego_names,ix] = sort(stego_names);
S = stego.F(ismember(stego.names,names),:);
S = S(ix,:);
% 至此,第i行stego matrix 与 低i行cover matrix是一对

%% 分割为training_set & testing_set
%
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1));
random_permutation = randperm(size(C,1));
if ~exist('numOfTrain', 'var')
  numOfTrain = round(size(C,1)/2);
end
training_set = random_permutation(1:numOfTrain);
testing_set = random_permutation(numOfTrain+1 : end);
if(numOfTrain>size(C,1)/2)
  testing_set = random_permutation;
end
training_names = names(training_set);
testing_names = names(testing_set);
% Prepare training features
TRN_cover = C(training_set,:);
TRN_stego = S(training_set,:);
% Prepare testing features
TST_cover = C(testing_set,:);
TST_stego = S(testing_set,:);

%% Train ensemble with all settings default - automatic search for the
% optimal subspace dimensionality (d_sub), automatic stopping criterion for
% the number of base learners (L), both PRNG seeds (for subspaces and
% bootstrap samples) are initialized randomly.
% 'd_sub', 1600, 
settings = struct('verbose',2); % 针对SRM特征设置d_sub=1600
t0=datetime('now');
[trained_ensemble, results]= ensemble_training(TRN_cover,TRN_stego,settings);
fprintf('\n训练时间:'); disp(datetime('now')-t0);
save('trained_ensemble.mat','trained_ensemble');
save('results.mat','results');
a=1;
% Resulting trained classifier is a cell array containing individual base
% learners. Variable 'results' contains some additional information like
% optimal found parameters or the progress of the OOB error during the
% search.

% plot the results of the search for optimal subspace dimensionality
%{
figure(1);
clf;plot(results.search.d_sub,results.search.OOB,'.b');hold on;
plot(results.optimal_d_sub,results.optimal_OOB,'or','MarkerSize',8);
xlabel('Subspace dimensionality');ylabel('OOB error');
legend({'all attempted dimensions',sprintf('optimal dimension %i',results.optimal_d_sub)});
title('Search for the optimal subspace dimensionality');

% plot the OOB progress with the increasing number of base learners (at the
% optimal value of subspace dimensionality).

figure(2);
clf;plot(results.OOB_progress,'.-b');
xlabel('Number of base learners');ylabel('OOB error')
title('Progress of the OOB error estimate');
%}

% Time to test the performance of the classifier on the testing set that
% contains unseen samples. Even though OOB is an unbiased error estimate,
% we used it as a feedback for determining d_sub and L. Therefore, its
% value is an optimistic estimate of the real testing error.

%% Testing phase
%{
test_results_cover = ensemble_testing(TST_cover,trained_ensemble);
test_results_stego = ensemble_testing(TST_stego,trained_ensemble);

% Predictions: -1 stands for cover, +1 for stego
false_alarms = sum(test_results_cover.predictions~=-1);
missed_detections = sum(test_results_stego.predictions~=+1);
num_testing_samples = size(TST_cover,1)+size(TST_stego,1);
PFA = false_alarms / size(TST_cover, 1);
PMD = missed_detections / size(TST_stego, 1);
PE = (PFA+PMD) / 2;
fprintf('PE: %.4f  PFA: %.4f  PMD: %.4f\n',PE,PFA,PMD);
%}

%% 重复10次训练及测试
%{
% To speed-up the following example, we fix d_sub to 300
t0=datetime('now');
PMD= zeros(10,1); PFA=zeros(10,1);
d_sub = results.optimal_d_sub; % d_sub=1600;
settings = struct('d_sub',d_sub,'verbose',2);
for seed = 1:10
  RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));
  random_permutation = randperm(size(C,1));
  training_set = random_permutation(1:round(size(C,1)/2));
  testing_set = random_permutation(round(size(C,1)/2)+1:end);
  TRN_cover = C(training_set,:);
  TST_cover = C(testing_set,:);
  TRN_stego = S(training_set,:);
  TST_stego = S(testing_set,:);
  
  [trained_ensemble,~] = ensemble_training(TRN_cover,TRN_stego,settings);
  test_results_cover = ensemble_testing(TST_cover,trained_ensemble);
  test_results_stego = ensemble_testing(TST_stego,trained_ensemble);

  false_alarms = sum(test_results_cover.predictions~=-1);
  missed_detections = sum(test_results_stego.predictions~=+1);
  num_testing_samples = size(TST_cover,1)+size(TST_stego,1);
  PFA(seed) = false_alarms / size(TST_cover, 1);
  PMD(seed) = missed_detections / size(TST_stego, 1);

  fprintf('PMD %i: %.4f\n',seed, PMD(seed));
end
PE = (PMD+PFA)/2;
disp('#10次测试的平均值:'); 
fprintf('PFA: %.4f+/-%.4f\n', mean(PFA),std(PFA));
fprintf('PMD: %.4f+/-%.4f\n', mean(PMD),std(PMD));
fprintf('PE:  %.4f+/-%.4f\n', mean(PE),std(PE));
fprintf('\n重复10次耗时: '); disp(datetime('now')-t0);
%}

%% 设置settings
%{
% First, we can fix the random subspace dimensionality (d_sub) in order to
% avoid the expensive search. This can be useful for a fast research
% feedback.

settings = struct('d_sub',300);
[~,results] = ensemble_training(TRN_cover,TRN_stego,settings);

% The number of base learners (L) can be also fixed:
settings = struct('d_sub',300,'L',30);
[~,results] = ensemble_training(TRN_cover,TRN_stego,settings);

% Note: even for the fixed training data, resulting trained ensemble can
% have slightly different performance due to the stochastic components of
% bagging and random subspaces. The following loop executes ensemble
% training 10 times and outputs the average OOB error and its standard
% deviation.
settings = struct('d_sub',300,'L',30);
OOB = zeros(1,10);
for i=1:10
    [~,results] = ensemble_training(TRN_cover,TRN_stego,settings);
    OOB(i) = results.optimal_OOB;
end
fprintf('# -------------------------\n');
fprintf('Average OOB error = %.5f (+/- %.5f)\n',mean(OOB),std(OOB));

% Sometimes it may be useful to manually set the PRNG seeds for generating
% random subspaces and feature subsets for bagging. This will make the
% results fully reproducible. The following code, for example, should give
% the OOB error 0.1138:
settings = struct('d_sub',300,'L',30,'seed_subspaces',5,'seed_bootstrap',73);
[~,results] = ensemble_training(TRN_cover,TRN_stego,settings);

% You can fully suppress the ensemble output by setting 'verbose' to 0:
settings = struct('d_sub',300,'L',30,'verbose',0);
[~,results] = ensemble_training(TRN_cover,TRN_stego,settings);
fprintf('OOB error = %.4f\n',results.optimal_OOB);

% Alternatively, you can suppress everyting BUT the last line of the
% ensemble output (the one with the results) by setting 'verbose' to 2:
settings = struct('d_sub',300,'L',30,'verbose',2);
[~,results] = ensemble_training(TRN_cover,TRN_stego,settings);
%}