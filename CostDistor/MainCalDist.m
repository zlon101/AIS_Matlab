cRoot = 'E:\astego\StandExpers\covers\';
sRoot = 'E:\astego\StandExpers\CZL\';
D_CZL = 0;

% 遍历所有**格式文件
format='pgm';
cDirs = dir([cRoot,'*.',format]);
sDirs = dir([sRoot,'*.',format]);
num = length(cDirs);
D_CZL = zeros(num,1,'single');
count=0;  old=''; t0=datetime('now');
for i = 1:num
  % names{i}=coverDirs(i).name;
  cPath = [cRoot, cDirs(i).name];
  sPath = [sRoot, sDirs(i).name];

  % 计算stego 与 cover 之间的失真值
  D_CZL(i) =  calcuDist(single(imread(cPath)),single(imread(sPath)));
    
  % 打印
  count=count+1;
  msg=sprintf('- count: %3d/%d',count,num);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end

clear i cPath sPath num count old t0 msg cRoot sRoot...
  cDirs sDirs format;
% save('SUNWD_04.mat','SUNWD_04');
% fprintf('\n耗时: ');  disp(datetime('now')-t0);