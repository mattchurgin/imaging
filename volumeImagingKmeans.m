function [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(numKmeans,expectedNumberVolumes,nChannels,nSeries,fractionOfVariance,smoothingFrames,imageThreshold)
% volumeImagingKmeans finds clusters in a set of a volume images (or 2D
% image times series)
%
% numKmeans is the number of clusters to find
% expectedNumberVolumes is how many volume acquisitions are expected to be acquired
% for each odor delivery.  Sometimes acquisition fails, resulting in an
% incorrect number of volumes for a given delivery.  In this case the
% program automatically fills in zeros for that volume.
% nChannels is the number of image channels saved by scanimage
% nSeries is number of volume series to expect. important for
% pre-allocation
% fractionOfVariance designates how many principal components of data to
% keep for calculating kmeans
% smoothingFrames is number of frames to smooth (temporal smoothing)
% imageThreshold is the gray-scale intensity value above which a pixel must
% have on average to be considered for clustering
% Matt Churgin, August 2018
warning('off')

if nargin<2
    expectedNumberVolumes=13; % default number of volume time points
end

if nargin<3
    nChannels=2; % default number of image channels
end

if nargin<4
    nSeries=13; % default number of odor presentations
end

if nargin<5
    fractionOfVariance=0.75; % fraction of variance to keep for kmeans
end

if nargin<6
    smoothingFrames=1; % frames to smooth (no assumption about time per frame)
end

if nargin<7
    imageThreshold=10; % default pixel percentile threshold
end

folder1 = uigetdir;
home1 = folder1;
cd(home1)
currFolders = dir(folder1);
currFolders=currFolders(3:end);
greenImages=cell(1,nSeries);
redImages=cell(1,nSeries);
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
            
            if i==1
                unwrappedImages=zeros(size(greenUnwrapped,1)*nSeries,size(greenUnwrapped,2));
                unwrappedRedImages=zeros(size(redUnwrapped,1)*nSeries,size(redUnwrapped,2));
            end
           
            unwrappedImages(((size(greenUnwrapped,1)*(i-1)):(size(greenUnwrapped,1)*i)-1)+1,:)=greenUnwrapped;
            unwrappedRedImages(((size(redUnwrapped,1)*(i-1)):(size(redUnwrapped,1)*i)-1)+1,:)=redUnwrapped;
            
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
        red = smooth3(red,'gaussian',[3 3 1]);
        greenImages{i}=green;
        redImages{i}=red;
        
        %pre-allocate temporary unwrapped image
        tempunwrapped=zeros(size(green,3),size(green,1)*size(green,2));
        tempunwrappedred=zeros(size(red,3),size(red,1)*size(red,2));
        for k=1:size(green,3)
            tempunwrapped(k,:)=reshape(green(:,:,k),1,size(green,1)*size(green,2));
            tempunwrappedred(k,:)=reshape(red(:,:,k),1,size(red,1)*size(red,2));
        end
        if i==1
            unwrappedImages=[];
            unwrappedRedImages=[];
        end
        unwrappedImages=[unwrappedImages; tempunwrapped];
        unwrappedRedImages=[unwrappedRedImages; tempunwrappedred];
        imageSize=[size(green,1) size(green,2) 1];
        
        display(['loaded image series ' num2str(i) ' of ' num2str(length(currFolders))])
    end
end
display(['images loaded.  time elapsed: ' num2str(toc) ' seconds'])

if smoothingFrames>1
    % apply temporal smoothing to unwrappedImage for kmeans calculation
    unwrappedImages=conv2(ones(1,smoothingFrames)/smoothingFrames,1,unwrappedImages,'same');
    display(['applied temporal smoothing (' num2str(smoothingFrames) ' frames)'])
end


% use ind2sub to add pixel locations to the time series
% this allows us to add pixel location as a factor for training kmeans
% may need to multiple pixel location by a lambda factor to weight it more
% highly since it is only represented as 3 ppoints in a much longer time
% series
% 
% lambda=100;
% [Iind,Jind,Kind] = ind2sub(imageSize,1:size(unwrappedImages,2));
% % weight and normalize pixel locations
% Iind=lambda*(Iind-mean(Iind))./std(Iind);
% Jind=lambda*(Jind-mean(Jind))./std(Jind);
% Kind=lambda*(Kind-mean(Kind))./std(Kind);
% 
% display(['pixel locations calculated'])

% find mean and std of unwrapped images and get pixel threshold
meanUnwrapped=mean(unwrappedImages);
stdUnwrapped=std(unwrappedImages); % use std to feature scale
imageThresh=prctile(meanUnwrapped,imageThreshold);
lowpixels=meanUnwrapped<imageThresh;

unwrappedImagesMeanSubtracted=zeros(size(unwrappedImages,1),size(unwrappedImages,2));
unwrappedImagesMeanSubtractedNorm=zeros(size(unwrappedImages,1),size(unwrappedImages,2));
for j=1:size(unwrappedImages,1)
    unwrappedImagesMeanSubtracted(j,:)=unwrappedImages(j,:)-meanUnwrapped; % subtract mean only
    unwrappedImagesMeanSubtractedNorm(j,:)=unwrappedImagesMeanSubtracted(j,:)./stdUnwrapped;
end
display(['normalized unwrapped image'])

meanGreenImages = reshape(meanUnwrapped,[imageSize(1),imageSize(2),imageSize(3)]);

% get mean red image
meanRedUnwrapped=mean(unwrappedRedImages);
meanRedImages = reshape(meanRedUnwrapped,[imageSize(1),imageSize(2),imageSize(3)]);

% remove pixels below intensity threshold
unwrappedImagesMeanSubtracted(:,lowpixels)=NaN;
unwrappedImagesMeanSubtractedNorm(:,lowpixels)=NaN;
display(['low intensity pixels removed'])

% add pixel locations to unwrappedimage
%unwrappedImagesMeanSubtracted=[unwrappedImagesMeanSubtracted; Iind; Jind; Kind];
%display(['pixel locations added to unwrapped image'])

% apply PCA to reduce dimensionality
tic
display(['computing principal components'])
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(unwrappedImagesMeanSubtractedNorm');
for j=1:length(EXPLAINED)
    explainedVariance(j)=sum(EXPLAINED(1:j))/sum(EXPLAINED);
end

[temp inds]=find(explainedVariance>fractionOfVariance);
ncomponents=inds(1)-1; % number of PCs needed to retain fractionOfVariance
unwrappedReduced=zeros(size(unwrappedImagesMeanSubtracted,2),ncomponents);
unwrappedReduced=unwrappedImagesMeanSubtractedNorm'*COEFF(:,1:ncomponents);
%unwrappedReduced=unwrappedImagesMeanSubtracted'*COEFF(:,1:ncomponents);
display(['keeping ' num2str(ncomponents) ' components, which explain ' num2str(fractionOfVariance*100) '% of variance, for kmeans calculation'])
display(['principal components computed.  time elapsed = ' num2str(toc) ' seconds'])

% Part 2:  apply kmeans
tic
display(['beginning kmeans calculation'])
maxIteration=150;
numReplicates=10;
useParallel=0;
% if using parallel processing toolbox
if useParallel
    pool = parpool;                      % Invokes workers
    stream = RandStream('mlfg6331_64');  % Random number stream
    options = statset('UseParallel',1,'UseSubstreams',1,...
        'Streams',stream);
    
    [kmeansOut,C,sumd,D] = kmeans(unwrappedReduced,numKmeans,'Options',options,'MaxIter',maxIteration,'Display','iter','Distance','sqeuclidean','Replicates',numReplicates);
else % not using parallel processing toolbox
    [kmeansOut,C,sumd,D] = kmeans(unwrappedReduced,numKmeans,'MaxIter',maxIteration,'Display','iter','Distance','sqeuclidean');
end

% restore original image shape
kmeansOut = reshape(kmeansOut,[imageSize(1),imageSize(2),imageSize(3)]);

display(['finished calculating kmeans.  time elapsed: ' num2str(toc) ' seconds'])

% save data in current directory
save(['rawKmeans_' num2str(numKmeans) '_clusters_' num2str(fractionOfVariance) 'fractionOfVarianceKept.mat'],'folder1','kmeansOut','meanGreenImages','meanRedImages','imageIs2d','expectedNumberVolumes','C','sumd','explainedVariance','numKmeans','nChannels')