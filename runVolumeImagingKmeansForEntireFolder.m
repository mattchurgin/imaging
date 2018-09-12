%% run VolumeImagingkKmeans.m for all flies in a directory
% Matt Churgin, September 2018
folderToRunKmeansOn='Volumes';

startDir=uigetfolder(); % user selects the folder containing each fly's image folders
currFolders = dir(startDir);
currFolders=currFolders(3:end);
for i=1:length(currFolders)
    actuallyAFolder(i)=currFolders(i).isdir;
end
currFolders(~actuallyAFolder)=[];

for i=1:length(currFolders)
    cd(currFolders(i).name)
    if exist('leftLobe')
        cd(['leftLobe\' folderToRunKmeansOn])
       [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(80,20,2,13,50,0.99);
        cd([startDir '\' currFolders(i).name])
    end
    if exist('rightLobe')
        cd(['rightLobe\'  folderToRunKmeansOn])
        [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(80,20,2,13,50,0.99);
        cd([startDir '\' currFolders(i).name])
    end
    cd(startDir)
end
disp(['done calculating kmeans.'])

%% run makeClusters.m and processVolumeImagingKmeans.m for all folders

folderToRunKmeansOn='Volumes';
rawKmeansOutput='rawKmeans_80clusters_0.99fractionOfVarianceKept_50replicates.mat';
clusterVolFile='clusterVols.mat';
volumeAcquisitionTime=1.2;
nChannels=2;

startDir=uigetfolder(); % user selects the folder containing each fly's image folders
currFolders = dir(startDir);
currFolders=currFolders(3:end);
for i=1:length(currFolders)
    actuallyAFolder(i)=currFolders(i).isdir;
end
currFolders(~actuallyAFolder)=[];

for i=1:length(currFolders)
    cd(currFolders(i).name)
    if exist('leftLobe')
        cd(['leftLobe\' folderToRunKmeansOn])
       [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
       [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
        cd([startDir '\' currFolders(i).name])
    end
    if exist('rightLobe')
        cd(['rightLobe\'  folderToRunKmeansOn])
       [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
       [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
        cd([startDir '\' currFolders(i).name])
    end
    cd(startDir)
end
disp(['done creating clusters and responses.'])