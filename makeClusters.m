function [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput)
% makeClusters takes the rawKmeansOutput file and creates clusterVolU and
% clusterInfoU cell arrays
% Matt Churgin, September 2018

load(rawKmeansOutput)

warning('off')

if numReplicates>1
    disp('processing multiple replicates')
    [clusterVolU clusterInfoU clusterArea nClustersFound] = processMultipleKmeans(rawKmeansOutput);
    disp('processed clusters for multiple replicates')
    [clusterVolConsensus clusterInfoConsensus] = calculateConsensusClustersOverlap(clusterVolU,clusterInfoU);
    disp('found consensus clusters')
    
    clear clusterVolU clusterInfoU
    clusterVolU=clusterVolConsensus;
    clusterInfoU=clusterInfoConsensus;
else
    disp('processing single replicate')

    [clusterVols clusterInfo clusterVolU clusterInfoU]=visualizeClusters(kmeansOut,400,20000,50);
end

showClusters(clusterVolU,clusterInfoU)

numClusters=length(clusterVolU);
save('clusterVolFile.mat','clusterVolU','clusterInfoU')