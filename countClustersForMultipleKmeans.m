%%
% run visualize clusters for multiple kmeans runs with multiple replicates
clear all
close all
load(uigetfile());
numKmeans=numKmeans(1:38);
nClustersFound=zeros(length(numKmeans),length(kmeansOut{1}));
clusterVolU=cell(length(numKmeans),length(kmeansOut{1}));
clusterInfoU=cell(length(numKmeans),length(kmeansOut{1}));
clusterArea=cell(length(numKmeans),length(kmeansOut{1}));
for j=1:length(numKmeans)
    for jj=1:length(kmeansOut{j})
    [clusterVol clusterInfo uniqueCV uniqueCI uniqueCI2d] = visualizeClusters((kmeansOut{j}{jj}),300,50000,50,0);
    clusterVolU{j}{jj}=uniqueCV;
    clusterInfoU{j}{jj}=uniqueCI;
    clusterInfo2dU{j}{jj}=uniqueCI2d;
    clusterArea{j}{jj}=zeros(1,length(uniqueCI));
    for z=1:length(uniqueCI)
       clusterArea{j}{jj}(z)=uniqueCI{z}.Area; 
    end
    nClustersFound(j,jj)=length(uniqueCV);
    disp(['found clusters for kmeans = ' num2str(numKmeans(j))])
    end
end
%%
% plot num clusters vs numkmeans
close all
figure

for j=1:length(numKmeans)
    for jj=1:size(nClustersFound,2)
   hold on
   plot(numKmeans(j),nClustersFound(j,jj),'ko')
    end
end
%errorbar(numKmeans,mean(nClustersFound,2),std(nClustersFound'),'LineWidth',2)
plot(numKmeans,mean(nClustersFound,2),'LineWidth',3)
xlabel('# k-means clusters')
ylabel('Final # clusters')
box off
set(gca,'FontSize',15)

mCA=zeros(size(nClustersFound,1),size(nClustersFound,2));
for j=1:size(nClustersFound,1)
    for jj=1:size(nClustersFound,2)
        mCA(j,jj)=mean(clusterArea{j}{jj});
    end
end

figure

for j=1:length(numKmeans)
    for jj=1:size(nClustersFound,2)
   hold on
   plot(numKmeans(j),mCA(j,jj),'ko')
    end
end
%errorbar(numKmeans,mean(nClustersFound,2),std(nClustersFound'),'LineWidth',2)
plot(numKmeans,mean(mCA,2),'LineWidth',3)
xlabel('# k-means clusters')
ylabel('Average Cluster Volume (pixels)')
box off
set(gca,'FontSize',15)