function PE=learnerTest(TST_cover,TST_stego,learner)
% ≤‚ ‘ learner
test_results_cover = ensemble_testing(TST_cover,{learner});
test_results_stego = ensemble_testing(TST_stego,{learner});

% Predictions: -1 stands for cover, +1 for stego
false_alarms = sum(test_results_cover.predictions~=-1);
missed_detections = sum(test_results_stego.predictions~=+1);
num_testing_samples = size(TST_cover,1)+size(TST_stego,1);
%% testing_error
PFA = false_alarms / size(TST_cover, 1);
PMD = missed_detections / size(TST_stego, 1);
PE = (PFA+PMD) / 2;
PE = round(PE,3); PFA=round(PFA,3); PMD=round(PMD,3);
fprintf('Testing error: %.4f\n',PE);