function [clusterVolU clusterInfoU grnResponse grnResponseNorm t] = processVolumeImagingKmeans(rawKmeansOutput,volumeAcquisitionTime,nChannels)
% processVolumeImagingKmeans takes the output of volumeImagingKmeans and
% prunes trivial clusters (removes clusters of too small or too large size)
% and produces the time series response of each cluster to each odor
% delivered
%
% rawKmeansOutput is a .mat file containing raw kmeans output
% volumeAcquisitionTime is the repetition rate of volume acquisition in
% seconds
% nChannels is the number of image channels saved by scanimage
% Matt Churgin, August 2018

tic
load(rawKmeansOutput)


% load images
home1 = pwd;
cd(home1)
currFolders = dir(home1);
currFolders=currFolders(3:end);
%greenImages=cell(1,nSeries);
%redImages=cell(1,nSeries);
imageIs2d=0; % track whether volume image or 2d image series

for i=1:length(currFolders)
    actuallyAFolder(i)=currFolders(i).isdir;
end
currFolders(~actuallyAFolder)=[];

tic
% Part 1: read in images and unwrap for kmeans
for i = 1:length(currFolders)
    % if volume imaging a time series
    if expectedNumberVolumes>1
        try
            [green greenUnwrapped red redUnwrapped]=readVolumeImageSeries(currFolders(i).name,expectedNumberVolumes,nChannels);
            greenImages{i}=green;
            redImages{i}=red;
            imageSize=[size(green,1) size(green,2) size(green,3)];
            display(['loaded volume series ' num2str(i) ' of ' num2str(length(currFolders))]) 
        catch
            display(['folder ' currFolders(i).name ' could not be read'])
        end
    else % if expected volumes = 1, we are dealing with a 2D time series
        % processing is slightly different for 2d image series
        imageIs2d=1;
        possibleImages=dir(currFolders(i).name);
        imageFiles=[];
        for j=1:length(possibleImages)
            try
                if strfind([currFolders(i).name '/' possibleImages(j).name],'tif')
                    imageFiles=[imageFiles j];
                end
            catch
                if strfind([currFolders(i).name '\' possibleImages(j).name],'tif')
                    imageFiles=[imageFiles j];
                end
            end
        end
        
        try
            [green red]=readVolumeImage([currFolders(i).name '/' possibleImages(imageFiles).name],nChannels);
        catch
            [green red]=readVolumeImage([currFolders(i).name '\' possibleImages(imageFiles).name],nChannels);
        end
        green = smooth3(green,'gaussian',[3 3 1]);
        greenImages{i}=green;
        redImages{i}=red;
       
        imageSize=[size(green,1) size(green,2) 1];
        
        display(['loaded image series ' num2str(i) ' of ' num2str(length(currFolders))])
    end
end
display(['images loaded.  time elapsed: ' num2str(toc) ' seconds'])

warning('off')
% use visualizeClusters to identify clusters with likely true glomerular
% sizes (remove noise)
%clusterVols=visualizeClusters(kmeansOut);

[clusterVols clusterInfo clusterVolU clusterInfoU clusterInfo2dU]=visualizeClusters(kmeansOut,500,20000,15);
drawnow

% dummyi=1;
% for i=1:length(clusterVols)
%     if sum(sum(sum(clusterVols{i})))>0 
%         clusterVolsNonzero{dummyi}=clusterVols{i};
%         clusterInfoNonzero{dummyi}=clusterInfo{i};
%         tempC{dummyi}=C(i,:);
%         clusterSum(dummyi)=sumd(i);
%         dummyi=dummyi+1;
%     end
% end
% 
% % sort by within cluster kmeans sum
% [vs is]=sort(clusterSum);
% clear clusterVols clusterInfo
% for i=1:length(is)
%     clusterVols{i}=clusterVolsNonzero{is(i)};
%     clusterInfo{i}=clusterInfoNonzero{is(i)};
%     newC{i}=tempC{is(i)};
% end
numClusters=length(clusterVolU);
% 
% % merge clusters with too similar centroids?
% % find distance between each cluster centroid
% for i=1:numClusters
%     for j=1:numClusters
%         cij(i,j)=sqrt(sum((newC{i}-newC{j}).^2));
%     end
% end


% [GloSig] = gloSig(ai_data,output,mask);
% [GloSig2] = gloSig(ai_data,output2,mask);

grnResponse=zeros(length(greenImages),numClusters,size(greenImages{1},4));
%redResponse=zeros(length(greenImages),numClusters,size(greenImages{1},4));
for i=1:length(greenImages)
    currBaseLine=nanmean(greenImages{i},4);
    for j = 1:numClusters
        for k=1:size(greenImages{1},4)
            % calculate df/f 
            grnResponse(i,j,k)=100*nanmean(nanmean(nanmean((greenImages{i}(:,:,:,k)-currBaseLine)./currBaseLine.*clusterVolU{j})));
            
            % calculate df/f using meanGreenImage
           %grnResponse(i,j,k)=100*mean(mean(mean((greenImages{i}(:,:,:,k)-meanGreenImages)./meanGreenImages.*(clusterVolU{j}))));
        end
    end
    display(['calculated cluster means for volumes ' num2str(i) ' of ' num2str(length(greenImages))])
end

grnResponseNorm=zeros(length(greenImages),numClusters,size(greenImages{1},4));
for i=1:length(greenImages)
    for j = 1:numClusters
        % subtract baseline
        grnResponseNorm(i,j,:)=grnResponse(i,j,:)-prctile(grnResponse(i,j,:),10);        
    end
    display(['calculated normalized response for volumes ' num2str(i) ' of ' num2str(length(greenImages))])
end

% calculate time vector
t=[1:size(greenImages{1},4)]*volumeAcquisitionTime;

figure
for i=1:size(grnResponseNorm,1)
    subplot(1,size(grnResponseNorm,1),i)
    
    imagesc(t,1:size(grnResponseNorm,2),squeeze(grnResponse(i,:,:)),[0 1])
    title(['odor ' num2str(i)])
    if i==1 
        xlabel('Time (s)')
    ylabel('Cluster #')
    title(['Air'])
    else 
        title(['Odor #' num2str(i-1)])
        set(gca,'YTick','')
    end
    set(gca,'FontSize',15)
end

% save data in current directory
save(['processedKmeans_' num2str(numKmeans) 'kmeans_' num2str(numClusters) 'uniqueclusters.mat'],'clusterInfoU','clusterVolU','clusterInfo2dU','grnResponse','grnResponseNorm','t')

disp(['finished. time to process raw k-means output: ' num2str(toc) ' seconds'])