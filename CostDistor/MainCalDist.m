coverRoot = 'E:\astego\Images\Experis\1covers\';
stegoRoot = 'E:\astego\Images\Experis\UNWD04\';
DUNWD = 0;

% 遍历所有**格式文件
coverDirs = dir([coverRoot, '*.bmp']);
stegoDirs = dir([stegoRoot, '*.bmp']);
num = length(coverDirs);
DUNWD = zeros(num,1,'single');
count=0;  old=''; t0=datetime('now');
for i = 1:num
  % names{i}=coverDirs(i).name;
  cPath = [coverRoot, coverDirs(i).name];
  sPath = [stegoRoot, stegoDirs(i).name];
  % 计算stego 与 cover 之间的失真值
  [DUNWD(i),R] =  calcuDist(single(imread(cPath)),single(imread(sPath)));
  %DHILL(i) = cacul_psnr(cPath, sPath);

  % 打印
  count=count+1;
  msg=sprintf('- count: %3d/%d',count,num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
clearvars -except DUNWD DHUGO DHILL;
% save('DHILL.mat','DHILL');
% fprintf('\n耗时: ');  disp(datetime('now')-t0);