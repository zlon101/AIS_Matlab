function analyze_histogram(Fc, Fs, nameAlg)
% 部分CCPEV特征分析
% 输入像素矩阵或者文件名
if(~exist('nameAlg', 'var'))
    nameAlg='';
end

fig=figure;
plot(Fc, 'k-');hold on;
plot(Fs, 'rx')
% title(nameAlg);
xlabel('CCPEV特征维度' ,'FontSize',14);
ylabel('特征值','FontSize',14);
legend({'cover','stego'}, 'Interpreter', 'none', 'FontSize',14)
ax = gca;
spreadAxes(ax);

ax.XLim = [1, length(Fc)];
xaxis = [1,11,66,165,193];
ax.XTick = [1:50:length(Fc), length(Fc)];
ax.YLim = [0,1];

fig.WindowStyle='normal';
fig.PaperPositionMode = 'manual';
fig.Units =  'pixels';
% statistic(coverQDCT, nameAlg, stegQDCT);
end

function F = subCcpev(QDCT)
% 提取部分CCPEV特征
% [F]=ExtractFeatures(image_name);

%----------------GLOBAL HISTOGRAM---------------------
% 值为-6:6的频数
H = hist(QDCT(:),-6:6);
F = H(2:end-1)/sum(H);
%--------------LOCAL HISTOGRAMs-----------------------
modes=[2 1;3 1;1 2;2 2;1 3];
for i=1:5
    Original = QDCT(modes(i,1):8:end,modes(i,2):8:end);
    H = hist(Original(:),-6:6);
    F = [F H(2:end-1)/sum(H)];
end
%---------------DUAL HISTOGRAM------------------------
%Dual histograms in range -5 to 5
modes=[2 1;1 2;3 1;2 2;1 3;4 1;3 2;2 3;1 4];
for value=-5:5
  T = zeros(1,9);
  for i=1:9
      T(i)=sum(sum(QDCT(modes(i,1):8:end, modes(i,2):8:end)==value));
  end
  F = [F T/max( sum(QDCT(:)==value), 1)];
end
end