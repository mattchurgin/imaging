function [greenChannel greenChannelUnwrapped redChannel redChannelUnwrapped] = readVolumeImageSeries(directory,expectedNumberVolumes,nChannels)
% readVolumeImageSeres reads a series of volume images acquired on two-photon microscope
% saves the green channel and red channel in separate image stacks
% Assumes all images in a directory are part of the time series and sorted
% Assumes .tif files
% nChannels is the number of image channels saved by scanimage
% Dependencies:  Uses readVolumeImage.m to read each volume

currfiles=dir(directory);
imageFiles=[];
for i=1:length(currfiles)
    if strfind(currfiles(i).name,'tif')
        imageFiles=[imageFiles i];
    end
end

if nargin<2
    numVolumes=length(imageFiles);
elseif length(imageFiles)==expectedNumberVolumes
    numVolumes=expectedNumberVolumes;
elseif length(imageFiles)<expectedNumberVolumes
    numVolumes=expectedNumberVolumes;
end

% get dimensions of each volume
try
    [gc1 rc1]=readVolumeImage([directory '/' currfiles(imageFiles(1)).name],nChannels);
catch
    [gc1 rc1]=readVolumeImage([directory '\' currfiles(imageFiles(1)).name],nChannels);
end

% parameters to smooth images
% eventually want to use size of the image to determine the smoothing parameters (higher pixel count means you can smooth more)
[xs ys zs]=size(gc1);
% default parameters
useDefaults=0;
if useDefaults
    gsSize=0.65;
    smoothingFilter=[3 3 3];
else
    gsSize=0.65;
    smoothingFilter=[5 7 3];
end

if gsSize>0
    gc1=smooth3(gc1,'gaussian',smoothingFilter,gsSize);
    rc1=smooth3(rc1,'gaussian',smoothingFilter,gsSize);
end

greenChannel=zeros(size(gc1,1),size(gc1,2),size(gc1,3),numVolumes);
redChannel=zeros(size(gc1,1),size(gc1,2),size(gc1,3),numVolumes);
greenChannelUnwrapped=zeros(numVolumes,size(gc1,1)*size(gc1,2)*size(gc1,3));
redChannelUnwrapped=zeros(numVolumes,size(gc1,1)*size(gc1,2)*size(gc1,3));

greenChannel(:,:,:,1)=gc1;
redChannel(:,:,:,1)=rc1;
greenChannelUnwrapped(1,:)=gc1(:);
redChannelUnwrapped(1,:)=rc1(:);

for i=2:numVolumes
    try
        try
            [gc1 rc1]=readVolumeImage([directory '/' currfiles(imageFiles(i)).name],nChannels);
        catch
            [gc1 rc1]=readVolumeImage([directory '\' currfiles(imageFiles(i)).name],nChannels);
        end
        
        if gsSize>0
            gc1=smooth3(gc1,'gaussian',smoothingFilter,gsSize);
            rc1=smooth3(rc1,'gaussian',smoothingFilter,gsSize);
        end
        
        greenChannel(:,:,:,i)=gc1;
        redChannel(:,:,:,i)=rc1;
        greenChannelUnwrapped(i,:)=gc1(:);
        redChannelUnwrapped(i,:)=rc1(:);
    catch
        try
            display(['volume ' num2str(i) ' could not be read in folder: ' directory ', filename: ' currfiles(imageFiles(i)).name '.  File likely corrupted.'])
        catch
            display(['volume ' num2str(i) ' could not be read in folder: ' directory ', filename does not exist'])
        end
    end
end