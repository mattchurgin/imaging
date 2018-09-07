% show clusters for random replicates of kmeans and the consensus
clear all
close all

% load raw kmeans output
load(uigetfile());
ntoshow=5;
temp=randperm(length(kmeansOut));
indsToShow=temp(1:ntoshow);

nClustersFound=zeros(length(numKmeans),length(kmeansOut{1}));
clusterVolU=cell(1,length(kmeansOut));
clusterInfoU=cell(1,length(kmeansOut));
clusterArea=cell(1,length(kmeansOut));

for jj=1:length(indsToShow)
    [clusterVol clusterInfo uniqueCV uniqueCI uniqueCI2d] = visualizeClusters((kmeansOut{indsToShow(jj)}),300,50000,50);

    showClusters(uniqueCV,uniqueCI)
    title(['Iteration #' num2str(indsToShow(jj))])
end

% load consensus data
load(uigetfile());
showClusters(clusterVolU,clusterInfoU)
title('Consensus')

%% calculate consensus by running kmeans on 50 pixel cluster assignments
% show clusters for random replicates of kmeans and the consensus
clear all
close all

% load raw kmeans output
load(uigetfile());
load(uigetfile());
ntoshow=5;
temp=randperm(length(kmeansOut));
indsToShow=temp(1:ntoshow);

nClustersFound=zeros(length(numKmeans),length(kmeansOut{1}));
clusterVolU=cell(1,length(kmeansOut));
clusterInfoU=cell(1,length(kmeansOut));
clusterArea=cell(1,length(kmeansOut));

for jj=1:length(indsToShow)
    [clusterVol clusterInfo uniqueCV uniqueCI uniqueCI2d] = visualizeClusters((kmeansOut{indsToShow(jj)}),300,50000,50);

    showClusters(uniqueCV,uniqueCI)
    title(['Iteration #' num2str(indsToShow(jj))])
end

% load consensus data

[clusterVol clusterInfo uniqueCV uniqueCI uniqueCI2d] = visualizeClusters((kmeansR),300,50000,50);
showClusters(uniqueCV,uniqueCI)
title('Consensus')