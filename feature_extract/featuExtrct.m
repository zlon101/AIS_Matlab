% 提取目录中所有图像的特征
%%
cDir= 'E:\astego\Images\standard_images\covers\';
sDir = 'E:\astego\Images\StandExpers\czl4\';
name= '195.pgm';
cPath=[sDir,name];
% S_CZL4_SRM_04 = getFeatures(sDir, -1);

%% 加载从C++程序中提取出来的特征数据
clear;
path1 = 'F:\astego\特征CPP\stegos\锐化_Am0.8_HUGO_03\';
% path2 = 'E:\astego\特征CPP\锐化_Am0.8\2\';
% path3 = 'E:\astego\特征CPP\HUGOOPT1_04\3\';
% path4 = 'E:\astego\特征CPP\HUGOOPT1_04\4\';

F1 = LoadFeature(path1);
% F2 = LoadFeature(path2);
% F3 = LoadFeature(path3);
% F4 = LoadFeature(path4);

S_SharpAm08_HUGO_03SRM.names= F1.names;
S_SharpAm08_HUGO_03SRM.F= single(F1.F);
clear F1 F2 F3 F4;
[S_SharpAm08_HUGO_03SRM.names, ind]= sort(S_SharpAm08_HUGO_03SRM.names);
S_SharpAm08_HUGO_03SRM.F = S_SharpAm08_HUGO_03SRM.F(ind, :);
save('S_SharpAm08_HUGO_03SRM', 'S_SharpAm08_HUGO_03SRM');
%}