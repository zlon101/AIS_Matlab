function effectDims = featEffect(Fc, Fs, rate)
% 分析diff=Fs-Fc,哪个分量差异最大,即哪个特征分量对隐写修改最敏感
% diff可以是单个样本，也可以是多个样本
% 样本中偏差最大的那个分量及其出现的次数
if(~exist('rate','var'))
    rate = 0.1;
end
nd = size(Fc,2);            % 维度
nMax = round(rate*nd);          % 最敏感的n个分量
diff = (Fs-Fc)./(Fc+eps);   % diff是矩阵
D = abs(diff);
[~,dims] = sort(D,2,'descend');    % 按行降序排列
components = dims(:, 1:nMax);
components = components(:);
[freqs,dims] = histFreqs(components);
[~,ind] = sort(freqs, 'descend');
dimSored = dims(ind);
effectDims = dimSored(1:nMax);
% saveImgDiffeature(diff(1:4), dirSavaImg, 'CCPEV');

end

function saveImgDiffeature(D, dirSavaImg, featureName)
% 保存每个样本的(Fs-Fc)./Fc
% 一行代表一个样本

numSamp = size(D,1);
count=0;
old='';
while(count<=numSamp)
    close all;fh=figure;    
    set(fh,'visible','off');
    for i=1:1
        count=count+1;
        if(count>numSamp)
            break;
        end
        %subplot(2,2,i);
        plot( D(count,:), '-kx');
        title([featureName,'-',int2str(count)]);
        xlabel('特征维度');
        ylabel('Fs-Fc');
    end  
    filename=[int2str(count),'.jpg'];
    saveas(fh,[dirSavaImg, filename]);
    msg=sprintf('- count: %3d/%d',count,numSamp);
    fprintf([repmat('\b',1,length(old)),msg]);
    old=msg;
end

end

function [freqs,ind]=histFreqs(D)
% 计算向量中频数最大的值和索引 freqs、ind
% D:整数,
% [freqs,edges,bin] = histcounts(D,'BinMethod','integers');
% ind = ceil( edges(1:end-1) );
D=D(:);
[freqs, ind]=hist(D, min(D):max(D));

end
