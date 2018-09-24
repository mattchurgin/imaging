%% run VolumeImagingkKmeans.m for all flies in a directory
% Matt Churgin, September 2018
clear all
folderToRunKmeansOn{1}='Volumes2';
%folderToRunKmeansOn{2}='Volumes2';
daysToProcess{1}='180831_pairedbehaviorandimaging';
daysToProcess{2}='180906';
daysToProcess{3}='180911_pairedbehaviorimaging';
daysToProcess{4}='180914_pairedbehaviorimaging';
fractionOfVarianceTokeep=0.75;
replicates=50;
nkmeans=80;
volsPerOdor=20;
nchannels=2;
nvolseries=13;

homeDir=pwd;
for days=1:length(daysToProcess)
    cd(daysToProcess{days})
    startDir=pwd;
    currFolders = dir(startDir);
    currFolders=currFolders(3:end);
    for i=1:length(currFolders)
        actuallyAFolder(i)=currFolders(i).isdir;
    end
    currFolders(~actuallyAFolder)=[];
    
    for i=1:length(currFolders)
        cd(currFolders(i).name)
        if exist('leftLobe')
            try
                for j=1:length(folderToRunKmeansOn)
                    cd(['leftLobe\' folderToRunKmeansOn{j}])
                    [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(nkmeans,volsPerOdor,nchannels,nvolseries,replicates,fractionOfVarianceTokeep);
                    cd([startDir '\' currFolders(i).name])
                end
            catch
            end
        end
        if exist('rightLobe')
            try
                for j=1:length(folderToRunKmeansOn)
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                    [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(nkmeans,volsPerOdor,nchannels,nvolseries,replicates,fractionOfVarianceTokeep);
                    cd([startDir '\' currFolders(i).name])
                end
            catch
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end
disp(['done calculating kmeans.'])

%% run makeClusters.m and processVolumeImagingKmeans.m for all folders

folderToRunKmeansOn{1}='Volumes2';
%folderToRunKmeansOn{2}='Volumes2';
daysToProcess{1}='180831_pairedbehaviorandimaging';
daysToProcess{2}='180906';
daysToProcess{3}='180911_pairedbehaviorimaging';
daysToProcess{4}='180914_pairedbehaviorimaging';
rawKmeansOutput='rawKmeans_80clusters_0.75fractionOfVarianceKept_50replicates.mat';
clusterVolFile='clusterVols.mat';
volumeAcquisitionTime=1.2;
nChannels=2;

homeDir=pwd;
for days=1:length(daysToProcess)
    cd(daysToProcess{days})
    startDir=pwd;
    currFolders = dir(startDir);
    currFolders=currFolders(3:end);
    for i=1:length(currFolders)
        actuallyAFolder(i)=currFolders(i).isdir;
    end
    currFolders(~actuallyAFolder)=[];
    
    for i=1:length(currFolders)
        cd(currFolders(i).name)
        
        if exist('leftLobe')
            try
                for j=1:length(folderToRunKmeansOn)
                    cd(['leftLobe\' folderToRunKmeansOn{j}])
                    [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
                    [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
                    cd([startDir '\' currFolders(i).name])
                end
            catch
            end
        end
        if exist('rightLobe')
            try
                for j=1:length(folderToRunKmeansOn)
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                    [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
                    [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
                    cd([startDir '\' currFolders(i).name])
                end
            catch
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end
disp(['done creating clusters and responses.'])