%% run VolumeImagingkKmeans.m for all flies in a directory
% Matt Churgin, September 2018
clear all
folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';
fractionOfVarianceTokeep=0.75;
replicates=50;
nkmeans=80;
volsPerOdor=20;
nchannels=2;
nvolseries=13;

startDir=uigetdir(); % user selects the folder containing each fly's image folders
currFolders = dir(startDir);
currFolders=currFolders(3:end);
for i=1:length(currFolders)
    actuallyAFolder(i)=currFolders(i).isdir;
end
currFolders(~actuallyAFolder)=[];

for i=1:length(currFolders)
    cd(currFolders(i).name)
    if exist('leftLobe')
        for j=1:length(folderToRunKmeansOn)
            cd(['leftLobe\' folderToRunKmeansOn{j}])
            [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(nkmeans,volsPerOdor,nchannels,nvolseries,replicates,fractionOfVarianceTokeep);
            cd([startDir '\' currFolders(i).name])
        end
    end
    if exist('rightLobe')
        for j=1:length(folderToRunKmeansOn)
            cd(['rightLobe\'  folderToRunKmeansOn{j}])
            [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(nkmeans,volsPerOdor,nchannels,nvolseries,replicates,fractionOfVarianceTokeep);
            cd([startDir '\' currFolders(i).name])
        end
    end
    cd(startDir)
end
disp(['done calculating kmeans.'])

%% run makeClusters.m and processVolumeImagingKmeans.m for all folders

folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';
rawKmeansOutput='rawKmeans_80clusters_0.75fractionOfVarianceKept_50replicates.mat';
clusterVolFile='clusterVols.mat';
volumeAcquisitionTime=1.2;
nChannels=2;

startDir=uigetdir(); % user selects the folder containing each fly's image folders
currFolders = dir(startDir);
currFolders=currFolders(3:end);
for i=1:length(currFolders)
    actuallyAFolder(i)=currFolders(i).isdir;
end
currFolders(~actuallyAFolder)=[];

for i=1:length(currFolders)
    cd(currFolders(i).name)
    if exist('leftLobe')
        for j=1:length(folderToRunKmeansOn)
            cd(['leftLobe\' folderToRunKmeansOn{j}])
            [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
            [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
            cd([startDir '\' currFolders(i).name])
        end
    end
    if exist('rightLobe')
        for j=1:length(folderToRunKmeansOn)
            cd(['rightLobe\'  folderToRunKmeansOn{j}])
            [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
            [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
            cd([startDir '\' currFolders(i).name])
        end
    end
    cd(startDir)
end
disp(['done creating clusters and responses.'])