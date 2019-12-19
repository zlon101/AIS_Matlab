function embedInRoot(coverRoot,stegoRoot)
% 根据隐写算法对目录中的图像进行隐写
% coverRoot       图像目录
% stegRoot      输出载密图像目录
% numSample:    设置样本个数
% addpath(genpath('./featureFile'));      % 添加文件路径
t0 = datetime('now');
% coverRoot = 'E:\astego\Images\BOSS_ALL\';
% stegoRoot = 'E:\astego\Images\stegos\自适应隐写CZL1_04\';
format = 'pgm';
payLoad = single(0.4);
numSample = -1;

dirs  = dir([coverRoot,'*.',format]);
% imgFiles(1)=[];imgFiles(1)=[];
names=cell(length(dirs),1);           % 图像名,不含全部路径
if (~exist('numSample', 'var') || numSample<1)
  numSample=length(dirs);
end
nImages = length(dirs);
old='';
for i = 1:nImages
  names{i}=dirs(i).name;
  cPath=[coverRoot,dirs(i).name];
  % 嵌入算法
  
  stego=embedAlgCZL(single(imread(cPath)), payLoad);
  imwrite(uint8(stego), [stegoRoot,names{i}],format);
  %stego = S_UNIWARD(imread(cPath), payLoad);
  %stego = HUGO(single(imread(cPath)), payload, single([1,1]));
  %stego = HUGO_like(imread(cPath), payload);
  %stego = S_UNIWARD(imread([coverRoot,Names{i}]), payLoad);
  %stego = HILL([coverRoot,Names{i}], payLoad);
  % J_UNIWARD([coverRoot,Names{i}], [stegRoot,Names{i}], single(payLoad));
  % F5([coverRoot,Names{i}], payLoad, [stegRoot,Names{i}]);
  % nsf5_simulation([coverRoot,Names{i}], [stegRoot,Names{i}], payLoad,seed);
  % steg =LSBM(I, payLoad);imwrite(uint8(steg), [stegRoot, Names{i}]);

  % 打印
  msg=sprintf('- count: %3d/%d',i,numSample);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;    
  if i>=numSample
    break;
  end
end
fprintf('\nnumbel of img: %d\n', i);
fprintf('\n耗时: '); disp(datetime('now')-t0);
end