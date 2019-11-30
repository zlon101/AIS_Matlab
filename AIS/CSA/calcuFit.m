function [fits,Memory]=calcuFit(Abs,embedParas,Memory)
% 计算抗体适应度
% 
%%
srcData = single(imread(embedParas.srcPath));
num = length(Abs);
fits = zeros(num,1);
old='';
for i=1:num
    K = getKey(Abs(i,:));
    if(Memory{1,1}==K)
        vals = Memory{1,2};
        fits(i) = vals(Memory{1,1}==K);
    else
    % 锐化
    % [sharpedData, ~] = sharpen(srcData, Abs{i});
    sharpedData =  laplace(srcData, Abs(i,:));
    % 隐写
    sharpedStegoData = HUGO_like(uint8(sharpedData), embedParas.payLoad);
    % imwrite(sharpedData, sharpedPath, 'pgm');
    % imwrite(uint8(sharpedStegoData),sharpedStegoPath, 'pgm');
    %% 提取特征   
    %fetuStruct = getFeatures(sharpedPath);  Fc2 = fetuStruct.F;
    %fetuStruct = getFeatures(sharpedStegoPath);  Fs2 = fetuStruct.F;
    %fits(i) = norm(Fs2- Fc2);
    %fits(i) = cacul_psnr(sharpedPath, sharpedStegoPath);
    fits(i) =  calcuDist(sharpedData, sharpedStegoData);
    Memory = [Memory];
    % 打印
    msg=sprintf('- count: %3d/%d',i,num);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
    % if--end
    end
% for--end
end
%}

%  -------------------测试Castro----------
%{
f = '1 * x .* sin(4 * pi .* x) - 1 * y.* sin(4 * pi .* y + pi) + 1';
[x,y] = meshgrid(-1:0.05:1,-1:0.05:1); vxp = x; vyp = y;
vzp = eval(f);  % 函数值
Abs = cell2mat(Abs);
x = Abs(:,1);
y = Abs(:,2);
fits = eval(f);
imprime(1,vxp,vyp,vzp,x,y,fits,1,1);
% -------------------测试Castro----------
%}
% clearvars -except fits Memory;
end