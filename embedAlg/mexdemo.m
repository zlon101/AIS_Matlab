function mexdemo()
root = 'E:\astego\Images\BOSS_ALL\';
name = '1013.pgm';
payLoad = single(0.4);

cover = single(imread([root,name]));

%% าþะด
[rhoP1,rhoM1] = CostCZL(cover)