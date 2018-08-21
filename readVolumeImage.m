function [greenChannel redChannel] = readVolumeImage(fileName, nChannels)
% readVolumeImage reads a volume image (or 2d image series) acquired on two-photon microscope
% saves the green channel and red channel in separate image stacks
% nChannels is the number of image channels saved by scanimage
% Matt Churgin, August 2018
imageInfo=imfinfo(fileName);
nFrames=length(imageInfo);
ySize=imageInfo(1).Width;
xSize=imageInfo(1).Height;
greenChannel=zeros(xSize,ySize,nFrames/nChannels);
redChannel=zeros(xSize,ySize,nFrames/nChannels);

for i=nChannels:nChannels:nFrames
    % images are saved top to bottom (piezo scans from top to bottom), but we would rather display them from
    % bottom to top, so store images in reverse order
    redChannel(:,:,nFrames/nChannels+1-i/nChannels)=imread(fileName,'Index',i-(nChannels-1));
    greenChannel(:,:,nFrames/nChannels+1-i/nChannels)=imread(fileName,'Index',i-(nChannels-2));
end