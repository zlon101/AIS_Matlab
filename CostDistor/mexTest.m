% function [rhoP1,rhoM1] = mexTest
D= getSRMQ1FromSRM(S_HILL_04_SRM.F) - getSRMQ1FromSRM(C_BOSS_SRM.F);
DHILL= zeros(size(D,1),1);
for i=1:size(D,1)
  DHILL(i)=norm(D(i,:));
end
clear D i ans S_SUNWD_04_SRM S_HUGO_04_SRM S_HILL_04_SRM;
%%
%{
root = 'E:\astego\Images\Experis\';
name = '195.bmp';
cPath = [root,'1covers\',name];
sPath = [root,'HILL04\',name];
cover = single(imread(cPath));
[rhoP1,rhoM1] = CostHUGO(cover);

P = 1./rhoM1;
figure('name','HUGO');histogram(P);
figure('name','HUGO');imshow(P,[]);
clear root name cPath sPath cover;
% figure('name','µÈÈ«_-1¸ÅÂÊ');imshow(P,[]);
%}