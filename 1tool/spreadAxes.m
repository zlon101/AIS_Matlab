function spreadAxes(ax)
% 扩展坐标轴,使其铺满整个窗口
% ax:当前坐标轴
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
% ax.XTick=1:20:size(diff,2);
ax.TickLength = [0.001 0.2];
end