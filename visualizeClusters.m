function [clusterVol clusterInfo uniqueCV uniqueCI] = visualizeClusters(kmeansOut,pixelThresh,pixelThreshHigh,maxObjects)
% visualizeClusters takes the output of kmeans and smooths each cluster to
% remove noise
% Also plots the clusters for aiding identification
% pixelThresh is the minimum size of a cluster of pixels
% pixelThreshHigh is the maximum size of a cluster of pixels
% Choose the pixel thresholds such that very small and very large clusters
% of pixels are removed.  These tend to be noise. 
% maxObjects specifies the maximum number of objects that a cluster can
% have.  Typically clusters with very large number of objects are noise.


if nargin<3
    pixelThresh=100; % should depend on size of image
    pixelThreshHigh=10000; % should depend on size of image
end
if nargin<4
   maxObjects=15; 
end


figure
view(3);
axis tight
camlight
lighting gouraud
%lighting phong
hold on

nclusters=max(max(max(kmeansOut)));
clusterVol=cell(1,nclusters); % initialize cell of clusters
nIslands=1;
for i=1:nclusters
    temp=(kmeansOut==i);
    if pixelThresh>0
        origtemp=temp;
        CC = bwconncomp(temp); % find connected voxels within the cluster
        numPixels = cellfun(@numel,CC.PixelIdxList);
        
        % remove very small voxel clusters
        idx= find(numPixels<pixelThresh);
        for j=1:length(idx)
            temp(CC.PixelIdxList{idx(j)}) = 0;
        end
        
        % remove very large voxel clusters
        idx2= find(numPixels>pixelThreshHigh);
        for j=1:length(idx2)
            temp(CC.PixelIdxList{idx2(j)}) = 0;
        end
        
        
        %
%         %smooth the cluster map
%             temp=smooth3(temp,'gaussian',[3]);
%         %
%         %         % re-binarize
%             temp=temp>0;
%         %
%         %         % erode the cluster map
%             se = strel('sphere',1);
%             temp = imerode(temp, se);
%         % %
        
        % remove small and large objects again
        CC = bwconncomp(temp);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        
        % remove very small voxel clusters
        idx= find(numPixels<pixelThresh);
        for j=1:length(idx)
            temp(CC.PixelIdxList{idx(j)}) = 0;
        end
        
        % remove very large voxel clusters
        idx2= find(numPixels>pixelThreshHigh);
        for j=1:length(idx2)
            temp(CC.PixelIdxList{idx2(j)}) = 0;
        end

    end
    
    % get processed cluster info
    CC = bwconncomp(temp);
    clusterInfo{i}=regionprops(CC,'basic');
    
    if length(clusterInfo{i})<maxObjects
        % store cluster map
        clusterVol{i}=temp;
        
        % create separate cluster array composed only of single pixel islands
        for uniqueIsland=1:CC.NumObjects
            tempLoneIsland=zeros(1,size(temp(:),1));
            tempLoneIsland(CC.PixelIdxList{uniqueIsland})=1;
            uniqueCluster=reshape(tempLoneIsland,CC.ImageSize);
            CD = bwconncomp(uniqueCluster);
            uniqueCV{nIslands}=uniqueCluster;
            uniqueCI{nIslands}=regionprops(CD,'basic');
            
            p2=patch(isosurface(uniqueCluster),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
            %p2=reducepatch(p2,0.5);
            isonormals(uniqueCluster,p2)
            
            nIslands=nIslands+1;
        end
    else
        clusterVol{i}=[];
        clusterInfo{i}=[];
    end
   
    %drawnow
    %pause
    display(['processed k-means cluster ' num2str(i)])
end
disp(['finished.  found ' num2str(nIslands-1) ' isolated pixel islands.'])