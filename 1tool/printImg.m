% 保存图像
fig = gcf;
figure1 = figure('OuterPosition',[1 1 487.333333333333 270.666666666667]);

fig.WindowStyle='normal';
fig.PaperPositionMode = 'manual';
fig.Units =  'pixels';
fig.PaperUnits = 'inches';
% fig.PaperPosition = [0, 0, 300, 100];
fig.Renderer = 'painters';
imgRoot = 'E:\桌面文件\桌面\';
print([imgRoot, 'HUGO修改概率'],'-dpng', '-painters','-r0');