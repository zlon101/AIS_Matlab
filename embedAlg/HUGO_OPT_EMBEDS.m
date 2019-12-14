function HUGO_OPT_EMBEDS(coverRoot, stegoRoot,payload)
% 根据隐写算法对目录中的图像进行隐写
% coverRoot       图像目录
% stegRoot      输出载密图像目录
%%
t0 = datetime('now');
payload = single(str2double(payload));
dirs  = dir([coverRoot, '*.pgm']);
nImages = length(dirs);
old='';
for i = 1:nImages                % 遍历结构体就可以一一处理图片了
  cPath=[coverRoot,dirs(i).name];
  stego = HUGO(single(imread(cPath)), payload);
  imwrite(uint8(stego), [stegoRoot,dirs(i).name], 'pgm');
  
  % 打印
  msg=sprintf('- count: %3d/%d',i,nImages);
  fprintf([repmat('\b',1,length(old)),msg]);
  old=msg;    
end
fprintf('\n耗时: '); disp(datetime('now')-t0);
end