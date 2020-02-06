function mexdemo()
root = 'E:\astego\Images\BOSS_ALL\1\';
name = '1013.pgm';
payLoad = single(0.4);

cover = single(imread([root,name]));

%% าþะด
[rhoP1,rhoM1] = CostCZL_7(cover);
a=1;