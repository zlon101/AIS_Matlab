function mexdemo()
root = 'E:\astego\Images\';
name = '1013.pgm';
payLoad = single(0.4);
cover = single(imread([root,name]));

%% าþะด
[rhoP1,rhoM1] = CostCZL(cover);