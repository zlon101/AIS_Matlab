coverRoot = 'E:\astego\Images\BOSS_ALL\';
% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']); % coverDirs(1)=[];coverDirs(1)=[];
num = length(coverDirs);
save('coverDirs','coverDirs');
bestAbs = cell(num,2);
old='';
t0=datetime('now');
for i = 1:2
    load coverDirs;
    cPath = [coverRoot, coverDirs(i).name];
    bestAbs{i,1} = coverDirs(i).name;
    save('bestAbs','bestAbs'); clear bestAbs;
    clear coverDirs bestAbs bestFits TAbs;
    
    [bestFits,TAbs,~,~] = CSA(cPath);
    TAbs = cell2mat(TAbs);
    [vmin,~] = min(bestFits); 
    inds = (bestFits==vmin);
    Ab = min(TAbs(inds));
    
    load('bestAbs.mat'); bestAbs{i,2} = Ab;
    % 打印
    msg=sprintf('- count: %3d/%d',i,num);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
    %fprintf('\n耗时: '); disp(datetime('now')-t0);
    if(mod(i,100)==0)
        save('bestAbs', 'bestAbs');
    end
end
fprintf('\n耗时: '); disp(datetime('now')-t0);
save('bestAbs', 'bestAbs');

% figure('name','克隆选择');  title('克隆选择');
% plot(bestFits);
% hold on;
% plot(meanFits);

%{
clc;close all;
Root = 'E:\astego\Images\test\';
name = '195';
srcPath = [Root, name, '.pgm'];
srcStegoPath= [Root, name,'_Stego.pgm'];
sharpedPath = [Root,name,'_Sharped.pgm'];
sharpedStegoPath = [Root,name,'_SharpedStego.pgm'];
payLoad = single(0.4);

srcData = single(imread(srcPath));
stegoData = single(imread(srcStegoPath));
sharpedData = single(imread(sharpedPath));
sharpedStegoData = single(imread(sharpedStegoPath));

diff= sharpedData - srcData;
srcPSNR = cacul_psnr(sharpedPath, srcPath);
% 特征
fetuStruct = getFeatures(srcPath);  Fc = fetuStruct.F;
fetuStruct = getFeatures(srcStegoPath);  Fs = fetuStruct.F;
fetuStruct = getFeatures(sharpedPath);  Fc2 = fetuStruct.F;
fetuStruct = getFeatures(sharpedStegoPath);  Fs2 = fetuStruct.F;
% 欧氏距离
D1 = norm(Fs- Fc);
D2 = norm(Fs2- Fc2);
% imshow
figure('name','src');  imshow(uint8(srcData));
figure('name','sharped');  imshow(uint8(sharpedData));
figure('name','sharpedStego');  imshow(uint8(sharpedStegoData));
%}