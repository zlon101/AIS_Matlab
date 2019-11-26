% 计算cover 的代价和失真
coverRoot = 'E:\astego\Images\covers\锐化G3载体_Am1\';
stegoRoot = 'E:\astego\Images\stegos\HUGO\锐化G3_Am1_HUGO_04\';

% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.pgm']); % coverDirs(1)=[];coverDirs(1)=[];
stegoDirs = dir([stegoRoot, '*.pgm']); % stegoDirs(1)=[];stegoDirs(1)=[];
% names = cell(length(coverDirs),1);           % 图像名,不含全部路径
num = length(coverDirs);
D3 = zeros(num,1,'single');
count=0;  old=''; t0=datetime('now');
for i = 1:num
    % names{i}=coverDirs(i).name;
    cPath = [coverRoot, coverDirs(i).name];
    sPath = [stegoRoot, stegoDirs(i).name];
    % 计算stego 与 cover 之间的失真值
    % G1(i) = cacul_psnr(cPath, sPath);
    [D3(i),resid] =  calcuDist(cPath, sPath);
    
    
    % 打印
    count=count+1;
    msg=sprintf('- count: %3d/%d',count,num);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;    
end
% save('HILL_04','HILL_04');
% fprintf('\n耗时: ');  disp(datetime('now')-t0);