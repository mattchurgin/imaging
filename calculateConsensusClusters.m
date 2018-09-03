function [clusterVolConsensus clusterInfoConsensus] = calculateConsensusClusters(clusterVolU,clusterInfoU)

distThresh=3; % maximum distance allowed to call the clusters a match
fracIn=0.5; % fraction of iterations a cluster has to be in to be kept
pixelThresh=100;

% extract centroid for each cluster in each iteration
clusterCentroids=cell(1,length(clusterInfoU));
clusterArea=cell(1,length(clusterInfoU));
for i=1:length(clusterInfoU)
    for j=1:length(clusterInfoU{i})
        clusterCentroids{i}=[clusterCentroids{i} clusterInfoU{i}{j}.Centroid'];
        clusterArea{i}=[clusterArea{i} clusterInfoU{i}{j}.Area];
    end
end

% find euclidean distance between centroids
for ii=1:length(clusterCentroids)
    for jj=1:length(clusterCentroids)
        for i=1:length(clusterCentroids{ii})
            for j=1:length(clusterCentroids{jj})
                % euclidean distance matrix
                distM{ii}{jj}(i,j)=sqrt(sum((clusterCentroids{ii}(:,i)-clusterCentroids{jj}(:,j)).^2));
                
                % area distance (volume difference)
                distA{ii}{jj}(i,j)=abs(clusterArea{ii}(i)-clusterArea{jj}(j));
            end
        end
    end
end

for ii=1:length(clusterCentroids)
    for jj=1:length(clusterCentroids)
        [hm vm]=min(distM{ii}{jj},[],2);
        mind(ii,jj)=sum(hm<distThresh);
        toremove=hm>=distThresh;
        hm(toremove)=-1;
        vm(toremove)=-1;
        minDist{ii}(:,jj)=hm;
        minInd{ii}(:,jj)=vm;
    end
end

for i=1:length(clusterCentroids)
    tokeep=sum(minInd{i}>0,2)/length(clusterCentroids);
    nottokeep=tokeep<fracIn;
    minInd{i}(nottokeep,:)=-1;
end

% enumerate unique consensur clusters
indsToLeaveOut=cell(1,length(clusterCentroids));
numConsensusClusters=1;
for i=1:length(clusterCentroids)
    inds=find(minInd{i}(:,i)>0);
    inds=setdiff(inds,indsToLeaveOut{i});
    if length(inds)>0
        for j=1:length(inds)
            consensusClusters{numConsensusClusters}=minInd{i}(inds(j),:);
            numConsensusClusters=numConsensusClusters+1;
        end
    end
    % remove already used clusters from consideration
    for j=1:length(clusterCentroids)
        temp=minInd{i}(inds,j);
        temp=temp(find(temp>0));
        indsToLeaveOut{j}=[indsToLeaveOut{j} temp'];
    end
end

% find consensus pixels for each consensus cluster
clusterVolConsensus=cell(1,length(consensusClusters));
clusterInfoConsensus=cell(1,length(consensusClusters));
for i=1:length(consensusClusters)
    temp=zeros(size(clusterVolU{1}{1},1),size(clusterVolU{1}{1},2),size(clusterVolU{1}{1},3));
    currthresh=0;
    for j=1:length(consensusClusters{i})
        if consensusClusters{i}(j)>0
            temp=temp+clusterVolU{j}{consensusClusters{i}(j)};
            currthresh=currthresh+1;
        end
    end
    
    temp(temp<currthresh)=0;
    
    CC = bwconncomp(temp); % find connected voxels within the cluster
    numPixels = cellfun(@numel,CC.PixelIdxList);
    
    % remove all but maximum-sized cluster
    for j=2:length(CC.PixelIdxList)
        temp(CC.PixelIdxList{j}) = 0;
    end
    
    clusterVolConsensus{i}=temp;
    CD = bwconncomp(temp);
    clusterInfoConsensus{i}=regionprops(CD,'all');
end