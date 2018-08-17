function [clusterVols clusterInfo grnResponse grnResponseNorm t] = process2DImagingKmeans(rawKmeansOutput,frameRate)
% process2DImagingKmeans takes the output of volumeImagingKmeans 
% and produces the time series response of each cluster to each odor
% delivered
%
% rawKmeansOutput is a .mat file containing raw kmeans output
% framerate is the repetition rate of volume acquisition in
% seconds
% nChannels is the number of image channels saved by scanimage

filesep='/';
load(rawKmeansOutput)

% load images
home1 = folder1;
cd(home1)
currFolders = dir(folder1);
currFolders=currFolders(3:end);
unwrappedImages=[];
imageIs2d=0; % track whether volume image or 2d image series

dummyi=1;
% Part 1: read in images and unwrap for kmeans
for i = 1:length(currFolders)
    if currFolders(i).isdir
        % if volume imaging a time series
        if expectedNumberVolumes>1
            try
                [green greenUnwrapped red redUnwrapped]=readVolumeImageSeries(currFolders(i).name,expectedNumberVolumes,nChannels);
                greenImages{dummyi}=green;
                %size(green)
                redImages{dummyi}=red;
                imageSize=[size(green,1) size(green,2) size(green,3)];
                
                %[outputTemp,output2Temp,s11,tempKmean,grn{cnt,1},green,ana,red] = imageProcessor;
                
                %         rpi = length(tempKmean(1,:))/25;
                %         kmean_=cat(2,kmean_,tempKmean);
                %         %kmean_=cat(2,kmean_,tempKmean(:,:,rpi*6:rpi*18));
                %         output{cnt}= outputTemp;
                %         output2{cnt} = output2Temp;
                %         green2 = cat(3,green2,green);
                %         red2 = cat(3,red2,red);
                %         cnt = cnt+1;
                %         if isstruct(ana) & isempty(ai_data)
                %             ai_data = ana;
                %         end
                %         cd(home1)
                %
                display(['loaded volume series ' num2str(i) ' of ' num2str(length(currFolders))])
                dummyi=dummyi+1;
            catch
                display(['folder ' currFolders(i).name ' could not be read'])
            end
        else % if expected volumes = 1, we are dealing with a 2D time series
            % processing is slightly different for 2d image series
            imageIs2d=1;
            possibleImages=dir(currFolders(i).name);
            imageFiles=[];
            for j=1:length(possibleImages)
                if strfind([currFolders(i).name filesep possibleImages(j).name],'tif')
                    imageFiles=[imageFiles j];
                end
            end
            
            [green red]=readVolumeImage([currFolders(i).name filesep possibleImages(imageFiles).name],nChannels);
            
            green = smooth3(green,'gaussian',[3,3,1]);
            greenImages{dummyi}=green;
            %size(green)
            redImages{dummyi}=red;
            
            
            imageSize=[size(green,1) size(green,2) 1];
            dummyi=dummyi+1;
            display(['loaded image series ' num2str(i) ' of ' num2str(length(currFolders))])
        end
    end
end


nclusters=(max(max(kmeansOut)));
clusterVols=cell(1,nclusters); % initialize cell of clusters
for i=1:nclusters
    temp=(kmeansOut==i);
    
    % get final cluster info
    CC = bwconncomp(temp);
    clusterInfo{i}=regionprops(CC,'basic');
    
    % store cluster map
    clusterVols{i}=temp;

    display(['calculated cluster volume ' num2str(i)])
end

dummyi=1;
for i=1:length(clusterVols)
    if (sum(sum(clusterVols{i})))>0 
        clusterVolsNonzero{dummyi}=clusterVols{i};
        clusterInfoNonzero{dummyi}=clusterInfo{i};
        dummyi=dummyi+1;
    end
end
clusterVols=clusterVolsNonzero;
clusterInfo=clusterInfoNonzero;
numClusters=length(clusterVols);

grnResponse=zeros(length(greenImages),numClusters,size(greenImages{1},3));
%redResponse=zeros(length(greenImages),numClusters,size(greenImages{1},4));
for i=1:length(greenImages)    
    for j = 1:numClusters
        for k=1:size(greenImages{1},3)
            grnResponse(i,j,k)=(mean(mean((greenImages{i}(:,:,k)-meanGreenImages)./meanGreenImages.*(clusterVols{j}))));
            %redResponse(i,j,k)=mean(mean(mean(redImages{i}(:,:,:,k).*(clusterVols{j}))));
        end
    end
    display(['calculated cluster means for odor ' num2str(i) ' of ' num2str(length(greenImages))])
end

grnResponseNorm=zeros(length(greenImages),numClusters,size(greenImages{1},3));
for i=1:length(greenImages)
    for j = 1:numClusters
        %grnResponseNorm{i}(j,:)=grnResponse{i}(j,:)-mean(grnResponse{i}(j,:));
        %grnResponseNorm{i}(j,:)=grnResponse{i}(j,:)-grnResponse{i}(j,1);
        
        % normalize all responses (not sure if we should do this or not)
        %grnResponseNorm(i,j,:)=grnResponse(i,j,:)/sum(squeeze(grnResponse(i,j,:)));
        
        % subtract baseline
        grnResponseNorm(i,j,:)=grnResponse(i,j,:)-prctile(grnResponse(i,j,:),10);
        %grnResponseNorm(i,j,:)=grnResponse(i,j,:)-mean(grnResponse(i,j,:));
    end
    display(['calculated normalized response for odor ' num2str(i) ' of ' num2str(length(greenImages))])
end

% calculate time vector
t=[1:size(greenImages{1},3)]*frameRate;

figure
%ii=1:length(clusterVols);
for i=1:size(grnResponseNorm,1)
    subplot(1,size(grnResponseNorm,1),i)
    title(['Odor # ' num2str(i)])
    imagesc(squeeze(grnResponse(i,:,:)))
end

figure
for j=1:numClusters
plot(squeeze(grnResponseNorm(1,j,:)),'LineWidth',2)
hold on
end
legend

save(['processedKmeans_' num2str(numKmeans) 'kmeans_' num2str(numClusters) 'uniqueclusters.mat'],'clusterInfo','clusterVols','grnResponse','grnResponseNorm','t')