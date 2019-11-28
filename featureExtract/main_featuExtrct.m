% 提取目录中所有图像的特征
%%
clc;close all;
imgPath = 'E:\astego\Images\BOSS_SUB\1.pgm';
% F = getFeatures(imgPath, 1, 0);

% norm(T.F(1,:) - T.F(2,:));
% save('S_HUGO04_SRM', 'S_HUGO04_SRM');


%% 加载从C++程序中提取出来的特征数据
path1 = 'E:\astego\特征CPP\stego_HILL_04\1\';
path2 = 'E:\astego\特征CPP\stego_HILL_04\2\';
path3 = 'E:\astego\特征CPP\stego_HILL_04\3\';
path4 = 'E:\astego\特征CPP\stego_HILL_04\4\';
path5 = 'E:\astego\特征CPP\stego_HILL_04\5\';
path6 = 'E:\astego\特征CPP\stego_HILL_04\6\';

F1 = LoadFeature(path1);
F2 = LoadFeature(path2);
F3 = LoadFeature(path3);
F4 = LoadFeature(path4);
F5 = LoadFeature(path5);
F6 = LoadFeature(path6);

S_HILL_04_SRM.names= [F1.names;F2.names; F3.names;F4.names;F5.names;F6.names];
S_HILL_04_SRM.F= single([F1.F; F2.F; F3.F;F4.F;F5.F;F6.F]);
clear F1 F2 F3 F4;
[S_HILL_04_SRM.names, ind]= sort(S_HILL_04_SRM.names);
S_HILL_04_SRM.F = S_HILL_04_SRM.F(ind, :);
save('S_HILL_04_SRM', 'S_HILL_04_SRM');