function [clusterVolConsensus clusterInfoConsensus] = calculateConsensusClustersOverlap(clusterVolU,clusterInfoU)
% takes clusterVolU (single pixel island clusters) for multiple replicates
% using a single kmeans value.  requires pre-processing of rawKmeans output
% Matt Churgin, September 2018

overlapThresh=0.6;
fracIn=0.6; %.5 default % fraction of iterations a cluster has to be in to be kept
pixelThresh=0;
pixelConsensus=0.6; % fraction of iterations a pixel needs to be part of the cluster to be counted as consensus

% extract centroid for each cluster in each iteration
clusterCentroids=cell(1,length(clusterInfoU));
clusterLinearSize=cell(1,length(clusterInfoU));
for i=1:length(clusterInfoU)
    for j=1:length(clusterInfoU{i})
        clusterCentroids{i}=[clusterCentroids{i} clusterInfoU{i}{j}.Centroid'];
        clusterLinearSize{i}=[clusterLinearSize{i} nthroot(clusterInfoU{i}{j}.Area,3)];
    end
end
disp('loaded centroids')

% find euclidean distance between centroids
for ii=1:length(clusterCentroids)
    for jj=1:length(clusterCentroids)
        for i=1:size(clusterCentroids{ii},2)
            for j=1:size(clusterCentroids{jj},2)
                % euclidean distance matrix
                %distM{ii}{jj}(i,j)=sqrt(sum((clusterCentroids{ii}(:,i)-clusterCentroids{jj}(:,j)).^2));
                distM{ii}{jj}(i,j)=sqrt(sum((clusterCentroids{ii}(:,i)-clusterCentroids{jj}(:,j)).^2))/mean([clusterLinearSize{ii}(i) clusterLinearSize{jj}(j)]);
                distOverlap{ii}{jj}(i,j)=length(intersect(clusterInfoU{ii}{i}.PixelIdxList,clusterInfoU{jj}{j}.PixelIdxList))/max(length(clusterInfoU{ii}{i}.PixelIdxList),length(clusterInfoU{jj}{j}.PixelIdxList));
            end
        end
    end
    disp(['calculated distances for replicate ' num2str(ii) ' of ' num2str(length(clusterCentroids))])
end

for ii=1:length(clusterCentroids)
    for jj=1:length(clusterCentroids)
        [hm vm]=max(distOverlap{ii}{jj},[],2);
        mind(ii,jj)=sum(hm>overlapThresh);
        toremove=hm<=overlapThresh;
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

% enumerate unique consensus clusters
indsToLeaveOut=cell(1,length(clusterCentroids));
numConsensusClusters=1;
for i=1:length(clusterCentroids)
    inds=find(minInd{i}(:,i)>0);
    inds=setdiff(inds,indsToLeaveOut{i});
    if length(inds)>0
        for j=1:length(inds)
            consensusClusters{numConsensusClusters}=minInd{i}(inds(j),:);
            consensusClustersNum(numConsensusClusters)=sum(consensusClusters{numConsensusClusters}>0);
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
    
    temp(temp<(pixelConsensus*currthresh))=0;
    temp=logical(temp);
    CC = bwconncomp(temp); % find connected voxels within the cluster
    numPixels = cellfun(@numel,CC.PixelIdxList);
    
    csize=cellfun(@length, CC.PixelIdxList);
    [m mindex]=max(csize);
    aaa=1:length(CC.PixelIdxList);
    
    %%remove all but maximum-sized cluster
    for j=setxor(aaa,mindex)
        temp(CC.PixelIdxList{j}) = 0;
    end
    
    clusterVolConsensus{i}=temp;
    CD = bwconncomp(temp);
    clusterInfoConsensus{i}=regionprops(CD,'all');
end

todelete=[];
for i=1:length(clusterInfoConsensus)
    if length(clusterInfoConsensus{i})==0 || clusterInfoConsensus{i}.Area<pixelThresh
        todelete=[todelete i];
    end
end
clusterInfoConsensus(todelete)=[];
clusterVolConsensus(todelete)=[];
disp('found clusters. now merging non-unique clusters')

% merge clusters that are very close to each other
clusterConsensusCentroids=[];
clusterConsensusLinearSize=[];
for j=1:length(clusterInfoConsensus)
    clusterConsensusCentroids=[clusterConsensusCentroids clusterInfoConsensus{j}.Centroid'];
    clusterConsensusLinearSize=[clusterConsensusLinearSize nthroot(clusterInfoConsensus{j}.Area,3)];
end

prcOverlap=zeros(size(clusterConsensusCentroids,2),size(clusterConsensusCentroids,2));
for ii=1:size(clusterConsensusCentroids,2)
    for jj=1:size(clusterConsensusCentroids,2)
        distMConsensus(ii,jj)=sqrt(sum((clusterConsensusCentroids(:,ii)-clusterConsensusCentroids(:,jj)).^2))/mean([clusterConsensusLinearSize(ii) clusterConsensusLinearSize(jj)]);
    
        prcOverlap(ii,jj)=sum(clusterVolConsensus{ii}(:).*clusterVolConsensus{jj}(:))/mean([clusterInfoConsensus{ii}.Area clusterInfoConsensus{jj}.Area]);
    end
end

tomerge=prcOverlap>0.75;
indstomerge=find(sum(tomerge,1)>1);

if length(indstomerge>0)
    
    covered=setxor(indstomerge,1:length(clusterConsensusCentroids));
    clear merges
    mergesToMake=0;
    for i=1:length(covered)
        mergesToMake=mergesToMake+1;
        merges{mergesToMake}=covered(i);
    end
    
    for i=1:length(indstomerge)
        if ~any(ismember(covered,indstomerge(i)))
            temp=find(tomerge(:,indstomerge(i)));
            covered=[covered temp'];
            
            mergesToMake=mergesToMake+1;
            merges{mergesToMake}=temp;
        end
    end
    
    % finally merge the clusters
    for i=1:length(merges)
        if length(merges{i})==1
            clusterVolConsensusMerged{i}=clusterVolConsensus{merges{i}};
            clusterInfoConsensusMerged{i}=clusterInfoConsensus{merges{i}};
        else
            temp=zeros(size(clusterVolU{1}{1},1),size(clusterVolU{1}{1},2),size(clusterVolU{1}{1},3));
            currthresh=0;
            for j=1:length(merges{i})
                temp=temp+clusterVolConsensus{merges{i}(j)};
                currthresh=currthresh+1;
            end
            
            temp(temp<(pixelConsensus*currthresh))=0;
            temp=logical(temp);
            
            CC = bwconncomp(temp); % find connected voxels within the cluster
            numPixels = cellfun(@numel,CC.PixelIdxList);
            
            csize=cellfun(@length, CC.PixelIdxList);
            [m mindex]=max(csize);
            aaa=1:length(CC.PixelIdxList);
            
            %%remove all but maximum-sized cluster
            for j=setxor(aaa,mindex)
                temp(CC.PixelIdxList{j}) = 0;
            end
            
            clusterVolConsensusMerged{i}=temp;
            CD = bwconncomp(temp);
            clusterInfoConsensusMerged{i}=regionprops(CD,'all');
        end
    end
    
    clear clusterVolConsensus clusterInfoConsensus
    clusterVolConsensus=clusterVolConsensusMerged;
    clusterInfoConsensus=clusterInfoConsensusMerged;
end


todelete=[];
for i=1:length(clusterInfoConsensus)
    if length(clusterInfoConsensus{i})==0 || clusterInfoConsensus{i}.Area<pixelThresh
        todelete=[todelete i];
    end
end
clusterInfoConsensus(todelete)=[];
clusterVolConsensus(todelete)=[];