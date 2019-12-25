% 提取目录中所有图像的特征
%%
% cRoot= 'E:\astego\Images\BOSS_ALL\';
% name= '195.pgm';
% cPath=[cRoot,name];
% 
% cover= imread(cPath);
% SRM = SRM({cPath});

%% 加载从C++程序中提取出来的特征数据
path1 = 'E:\astego\标准图像集实验\特征CPP\SUNWD_04\';
% path2 = 'E:\astego\特征CPP\HUGOOPT1_04\2\';
% path3 = 'E:\astego\特征CPP\HUGOOPT1_04\3\';
% path4 = 'E:\astego\特征CPP\HUGOOPT1_04\4\';

F1 = LoadFeature(path1);
% F2 = LoadFeature(path2);
% F3 = LoadFeature(path3);
% F4 = LoadFeature(path4);

S_SUNWD_04_SRM.names= F1.names;
S_SUNWD_04_SRM.F= single(F1.F);
clear F1 F2 F3 F4;
[S_SUNWD_04_SRM.names, ind]= sort(S_SUNWD_04_SRM.names);
S_SUNWD_04_SRM.F = S_SUNWD_04_SRM.F(ind, :);
save('S_SUNWD_04_SRM', 'S_SUNWD_04_SRM');
%}