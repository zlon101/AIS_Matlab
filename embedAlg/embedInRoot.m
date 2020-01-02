function embedInRoot(coverRoot,stegoRoot,startInd,endInd)
% 根据隐写算法对目录中的图像进行隐写
% coverRoot       图像目录
% stegRoot      输出载密图像目录
% numSample:    设置样本个数

% coverRoot = 'E:\astego\StandExpers\covers\';
% stegoRoot = 'E:\astego\StandExpers\CZL\';
format = 'pgm';
payLoad = single(0.4);
dirs  = dir([coverRoot,'*.',format]);
nImages = length(dirs);
if(exist('startInd','var') && str2double(startInd)>0)
  startInd = single(str2double(startInd));
else
  startInd = 1;
end
if(exist('endInd','var'))
  endInd=single(str2double(endInd));
else 
  endInd=single(nImages);
end
% names=cell(length(dirs),1);
fprintf('# start\n#count: %d - %d\n',startInd,endInd);

old=''; t0 = datetime('now');
for i=startInd : endInd
  cPath=[coverRoot,dirs(i).name]; %names{i}=dirs(i).name;
  
  % 嵌入算法
  t0=tic;
  stego=embedAlgCZL(cPath, payLoad);
  disp(toc(t0));
  imwrite(uint8(stego), [stegoRoot,dirs(i).name],format);
  %stego = MiPOD( single(imread(cPath)), payLoad);
  %stego = HILL(cPath, payLoad);
  %stego = HUGO_like(imread(cPath), payLoad);
  %stego = HUGO(cPath, payLoad, []);
  %stego = S_UNIWARD(imread(cPath), payLoad);
  % J_UNIWARD([coverRoot,Names{i}], [stegRoot,Names{i}], single(payLoad));
  % steg =LSBM(I, payLoad);imwrite(uint8(steg), [stegRoot, Names{i}]);
  %stego = MG( single(imread(cPath)), payLoad );
  
  % 打印
  msg=sprintf('- count: %3d/%d',i,nImages);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;
end
fprintf('\nnumbel of img: %d\n', i);
fprintf('\n耗时: '); disp(datetime('now')-t0);
end