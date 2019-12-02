function embedInRoot(coverRoot, payLoad, stegoRoot, numSample)
% 根据隐写算法对目录中的图像进行隐写
% coverRoot       图像目录
% stegRoot      输出载密图像目录
% numSample:    设置样本个数
% addpath(genpath('./featureFile'));      % 添加文件路径
t0 = datetime('now');
coverRoot = 'E:\astego\Images\Experis\covers\';
stegoRoot = 'E:\astego\Images\Experis\stegos\UNWD\';
payLoad = single(0.3);
numSample = -1;

imgFiles  = dir([coverRoot, '*.pgm']);       % 遍历所有jpg格式文件
% imgFiles(1)=[];imgFiles(1)=[];
Names=cell(length(imgFiles),1);           % 图像名,不含全部路径
if (~exist('numSample', 'var') || numSample<1)
    numSample=length(imgFiles);
end
nImages = length(imgFiles);
old='';
for i = 1:nImages                % 遍历结构体就可以一一处理图片了
    imgName=[coverRoot imgFiles(i).name];
    Names{i}=imgFiles(i).name;
    % 嵌入算法
    stego = S_UNIWARD(imread([coverRoot,Names{i}]), payLoad);
    %stego = HILL([coverRoot,Names{i}], payLoad);
    imwrite(uint8(stego), [stegoRoot,Names{i}], 'pgm');
    % J_UNIWARD([coverRoot,Names{i}], [stegRoot,Names{i}], single(payLoad));
    % F5([coverRoot,Names{i}], payLoad, [stegRoot,Names{i}]);
    % nsf5_simulation([coverRoot,Names{i}], [stegRoot,Names{i}], payLoad,seed);
    % I = double(imread(imgName));
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