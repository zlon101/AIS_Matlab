% function [rhoP1,rhoM1] = mexTest
root = 'E:\astego\Images\Experis\';
name = '195.pgm';
cPath = [root,name];
sPath = [root,'HILL04\',name];
I = single(imread(cPath));

[rhoP1,rhoM1] = CostHUGO_like(I);
P = 1./rhoP1;
figure;imshow(P,[]);title('HUGOÐÞ¸Ä¸ÅÂÊ');
clear root name cPath sPath;

%{
D= getSRMQ1FromSRM(S_HILL_04_SRM.F) - getSRMQ1FromSRM(C_BOSS_SRM.F);
DHILL= zeros(size(D,1),1);
for i=1:size(D,1)
  DHILL(i)=norm(D(i,:));
end
%}