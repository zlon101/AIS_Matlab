% 计算cover 的代价和失真
coverRoot = 'E:\astego\Images\covers\锐化T1Cover\';
stegoRoot = 'E:\astego\Images\stegos\HUGO\锐化T1_HUGO04\';

% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']);
stegoDirs = dir([stegoRoot, '*.pgm']);
% names = cell(length(coverDirs),1); % 图像名,不含全部路径
num = length(coverDirs);
T1 = zeros(num,1,'single');
count=0;  old=''; t0=datetime('now');
for i = 1:num
    % names{i}=coverDirs(i).name;
    cPath = [coverRoot, coverDirs(i).name];
    sPath = [stegoRoot, stegoDirs(i).name];
    % 计算stego 与 cover 之间的失真值
    % G1(i) = cacul_psnr(cPath, sPath);
    [T1(i),resid] =  calcuDist(single(imread(cPath)),single(imread(sPath)));
    
    % 打印
    count=count+1;
    msg=sprintf('- count: %3d/%d',count,num);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;    
end
save('锐化T1.mat','T1');
% fprintf('\n耗时: ');  disp(datetime('now')-t0);