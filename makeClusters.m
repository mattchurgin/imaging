function [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput)
% makeClusters takes the rawKmeansOutput file and creates clusterVolU and
% clusterInfoU cell arrays
% Matt Churgin, September 2018

load(rawKmeansOutput)

topkmeans=5;

warning('off')

if numReplicates>1
    disp('processing multiple replicates')
    [clusterVolU clusterInfoU clusterArea nClustersFound] = processMultipleKmeans(rawKmeansOutput);
    disp('processed clusters for multiple replicates')
    
    % find best kmeans replicates and only use those for conenssus
    % calculation
    Csums = cellfun(@mean, sumd);
    [v i]=sort(Csums,'ascend');
    
    [clusterVolConsensus clusterInfoConsensus] = calculateConsensusClustersOverlap(clusterVolU(i(1:topkmeans)),clusterInfoU(i(1:topkmeans)));
    disp('found consensus clusters')
    
    clear clusterVolU clusterInfoU
    clusterVolU=clusterVolConsensus;
    clusterInfoU=clusterInfoConsensus;
else
    disp('processing single replicate')

    [clusterVols clusterInfo clusterVolU clusterInfoU]=visualizeClusters(kmeansOut,400,20000,50);
end

showClusters(clusterVolU,clusterInfoU)

save('clusterVols.mat','clusterVolU','clusterInfoU')