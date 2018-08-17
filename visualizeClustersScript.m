% script for 
clear all
close all

load(uigetfile) % load the output of volumeImagingKmeans
% convert kmeans output into individual volumes
figure
view(3);
axis tight
camlight
lighting gouraud
hold on
pixelThresh=100; % should depend on size of image
pixelThreshHigh=10000; % should depend on size of image
nclusters=max(max(max(kmeansOut)));
clusterVol=cell(1,nclusters); % initialize cell of clusters
for i=1:nclusters
    temp=(kmeansOut==i);
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
    
    % smooth the cluster map
    temp=smooth3(temp,'gaussian',3);
    
    % erode the cluster map
    se = strel('cube',3);
    temp = imerode(temp, se);
    
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
    
    % get final cluster info
    CC = bwconncomp(temp);
    clusterInfo{i}=regionprops(CC,'all');
    
    % store cluster map
    clusterVol{i}=temp;
    
    % visualize cluster
    p2=patch(isosurface(temp),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.5);
    isonormals(temp,p2)
    
    %drawnow
    %pause
    display(['calculated cluster volume ' num2str(i)])
end

