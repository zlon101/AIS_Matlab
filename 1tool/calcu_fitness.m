function [fit, diff, fits]=calcu_fitness(Fc, Fs)
% 适应度计算
%% 
numSamp = size(Fc,1);
% 单个样本特征
if(numSamp == 1)    
    [fit,diff] = singSample(Fc,Fs);
else    
    fits = zeros(numSamp, 1);
    for i=1:numSamp        
        fits(i) = singSample(Fc(i,:), Fs(i,:));
    end
    % 剔除异常样本
    ind= find(fits > 1e2);
    fits(ind)=0;    
    Fc(ind, :)=0;
    Fs(ind, :)=0;
    
    fit = mean(fits);    
    diff = mean(Fs-Fc, 1);    
end
end

function [fit, diff] = singSample(Fc, Fs)
%% 单个样本的适应度
%  方案1: 欧式距离
%  diff = ( Fs - Fc) ./ (Fc+eps);
%  ind = isinf(diff);
%  diff(ind)=Fs(ind);
diff = Fs-Fc;
% fitness = norm(diff);

%  方案2: FLD投影距离
load('learner.mat');    load('CMuSigma.mat');
C = bsxfun(@minus, Fc, CMuSigma.mu);  C = bsxfun(@rdivide, C, CMuSigma.sigma);
vC = C * learner.w;
S = bsxfun(@minus, Fs, CMuSigma.mu);  S = bsxfun(@rdivide, S, CMuSigma.sigma);
vS = S * learner.w;
fit = abs(vS - vC);
fit = round(fit, 3);        % 保留3为小数, 值越小越好
end