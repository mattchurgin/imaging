function [clusterVol clusterInfo uniqueCV uniqueCI uniqueCI2d] = visualizeClusters(kmeansOut,pixelThresh,pixelThreshHigh,maxObjects,plotYorN)
% visualizeClusters takes the output of kmeans and smooths each cluster to
% remove noise
% Also plots the clusters for aiding identification
% pixelThresh is the minimum size of a cluster of pixels
% pixelThreshHigh is the maximum size of a cluster of pixels
% Choose the pixel thresholds such that very small and very large clusters
% of pixels are removed.  These tend to be noise.
% maxObjects specifies the maximum number of objects that a cluster can
% have.  Typically clusters with very large number of objects are noise.
% plotYorN is a logical (0 or 1) designating whether to plot the clusters
% Matt Churgin, August 2018

if nargin<3
    pixelThresh=200; % should depend on size of image
    pixelThreshHigh=20000; % should depend on size of image
end
if nargin<4
    maxObjects=30;
end
if nargin<5
    plotYorN=1;
end

if plotYorN
    figure
    view(3);
    axis tight
    camlight
    lighting gouraud
    hold on
end

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
            uniqueClusterInfo=regionprops(CD,'basic');
            
            % look at z-centroid slice and get 2-d image properties to
            % determine if the cluster is real
            tempSliceCheck=bwconncomp(uniqueCluster(:,:,round(uniqueClusterInfo.Centroid(3))));
            if tempSliceCheck.NumObjects<100 % THERE IS A SMARTER WAY TO DO THIS (USE AREA OT DETERMINE NUMBER OF OBJECTS PERMISSIBLE????)
                
                uniqueCV{nIslands}=uniqueCluster;
                uniqueCI{nIslands}=uniqueClusterInfo;
                
                % save properties of 2d projected image
                
                CD2d = bwconncomp(sum(uniqueCluster,3)>0);
                nIslands
                uniqueClusterInfo2d=regionprops(CD2d,'all')
                uniqueCI2d{nIslands}=uniqueClusterInfo2d;
                
                if plotYorN
                    p2=patch(isosurface(uniqueCluster),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
                    isonormals(uniqueCluster,p2)
                    %
                                     %drawnow
                    %                 uniqueCI{nIslands}.Centroid
                                    % pause
                end
                nIslands=nIslands+1;
            end
        end
    else
        clusterVol{i}=[];
        clusterInfo{i}=[];
    end
    
    %display(['processed k-means cluster ' num2str(i)])
end
disp(['finished.  found ' num2str(nIslands-1) ' isolated pixel islands.'])