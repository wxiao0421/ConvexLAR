function []=plot_path(TP,FP,plotName)
% Function:
%   Plot solution paths
%
% Arguments:
%   FP: solution path at all grid points
%   TP: solution path at transition points 
%   plotName: title of the plot
% Output:
%   []

figure('Position',[0 0 700 400]);
title(plotName, 'FontWeight', 'bold', 'FontSize', 12);
p=(size(FP,1)-1)/2;
w=TP(1:p,:);
temp=sum(abs(w),1)/sum(abs(w(:,end)));
plotcol='rbgcm'
plotsty(1).sty='-';
plotsty(2).sty='-.'; 
plotsty(3).sty='--';
for i=1:size(w,1)
    line(sum(abs(FP(1:p,:)),1)/sum(abs(FP(1:p,end)),1), FP(i,:), 'Color', ...
        plotcol(mod(i,5)+1), 'LineStyle', plotsty(mod(i, 3)+1).sty)
end
t=axis;
hold on
for i=1:length(temp)
    line([temp(i) temp(i)], [t(3) t(4)], 'Color','k','LineWidth',1, ...
        'LineStyle', ':')
end
text(-0.10, 0, '\beta','FontSize',10, 'HorizontalAlignment', 'center',...
    'VerticalAlignment','middle', 'Rotation', 90)
xlabel('|\beta|/max|\beta|','FontSize', 10)
ax1 = gca;
set(ax1,'XColor','k','YColor','k', 'XAxisLocation','bottom','YAxisLocation',...
    'left','FontSize', 10, 'XLim', [0,1]);
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
set(ax2,'XTick',temp, 'XTickLabel',[0:1:1000] ,...
    'YTick',[],'FontSize', 10)
ytname=[1:1:1000];
[dump,ind]=sort(FP(1:p,end));
line([1 1], t(3:4), 'Color','k','Parent',ax2)
line([0 0], t(3:4), 'Color','k','Parent',ax2)
if sum(FP(ind,end)==0)==0
    set(ax2,'YTick',FP(ind,end), 'YTickLabel', ytname(ind),'FontSize', 10)
end