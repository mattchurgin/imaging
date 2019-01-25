clear all
close all

x=-10:0.05:10;
gs=2.5;
gau=1/sqrt(2*pi*gs^2)*exp(-(x.^2)/(2*gs^2));

figure
plot(x,gau,'k','LineWidth',3)
ylabel('probability')
xlabel('phenotype')
set(gca,'XTickLabel','\mu')
set(gca,'XTick',[0])
set(gca,'YTickLabel','')
set(gca,'YTick','')
box off
set(gca,'FontSize',20)