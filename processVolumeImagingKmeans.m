function [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels)
% processVolumeImagingKmeans takes the output of volumeImagingKmeans and
% prunes trivial clusters (removes clusters of too small or too large size)
% and produces the time series response of each cluster to each odor
% delivered
%
% clusterVolfile is a filename containing the clusterVolU and clusterInfoU
% cell arrays 
% rawKmeansOutput is the rawKmeansoutput file
% volumeAcquisitionTime is the repetition rate of volume acquisition in
% seconds
% nChannels is the number of image channels saved by scanimage
% Matt Churgin, August 2018
load(clusterVolFile)
load(rawKmeansOutput)

numClusters=length(clusterVolU);

tic    

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

grnResponse=zeros(length(greenImages),numClusters,size(greenImages{1},4));
%redResponse=zeros(length(greenImages),numClusters,size(greenImages{1},4));
for i=1:length(greenImages)
    currBaseLine=nanmean(greenImages{i},4);
    for j = 1:numClusters
        for k=1:size(greenImages{1},4)
            % calculate df/f 
            grnResponse(i,j,k)=100*nanmean(nanmean(nanmean((greenImages{i}(:,:,:,k)-currBaseLine)./currBaseLine.*(clusterVolU{j}))));
            
            % calculate df/f using meanGreenImage
           %grnResponse(i,j,k)=100*mean(mean(mean((greenImages{i}(:,:,:,k)-meanGreenImages)./meanGreenImages.*(clusterVolU{j}))));
        end
    end
    display(['calculated cluster means for volumes ' num2str(i) ' of ' num2str(length(greenImages))])
end

% calculate time vector
t=[1:size(greenImages{1},4)]*volumeAcquisitionTime;

figure
for i=1:size(grnResponse,1)
    subplot(1,size(grnResponse,1),i)
    
    imagesc(t,1:size(grnResponse,2),squeeze(grnResponse(i,:,:)),[0 1])
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
save(['clusterResponses.mat'],'clusterInfoU','clusterVolU','grnResponse','t')

disp('done')