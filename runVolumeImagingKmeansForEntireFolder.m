%% run VolumeImagingkKmeans.m for all flies in a directory
% Matt Churgin, September 2018
clear all

folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';

daysToProcess{1}='190211_orcoOX5gcamp_behaviorandimaging_copy';


fractionOfVarianceTokeep=.75;
ORNorPN=0; % 0 for processing ORNs, 1 for PNs
replicates=50;
if ORNorPN==1 % change to dynamically look for gh146 or orco in filename
    nkmeans=80; % use k=80 for processing PNs
else
    nkmeans=30; % use k=30 for processing ORNs
end
volsPerOdor=20;
nchannels=2;
nvolseries=13;

homeDir=pwd;
for days=1:length(daysToProcess)
    cd(daysToProcess{days})
    startDir=pwd;
    currFolders = dir(startDir);
    currFolders=currFolders(3:end);
    clear actuallyAFolder
    for i=1:length(currFolders)
        actuallyAFolder(i)=currFolders(i).isdir;
    end
    currFolders(~actuallyAFolder)=[];
    
    for i=1:length(currFolders)
        cd(currFolders(i).name)
        if exist('leftLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['leftLobe\' folderToRunKmeansOn{j}])
                    cd(['leftLobe\' folderToRunKmeansOn{j}])
                    [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(nkmeans,volsPerOdor,nchannels,nvolseries,replicates,fractionOfVarianceTokeep);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        if exist('rightLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['rightLobe\'  folderToRunKmeansOn{j}])
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                    [greenImages redImages kmeansOut unwrappedImages unwrappedReduced meanGreenImages meanRedImages] = volumeImagingKmeans(nkmeans,volsPerOdor,nchannels,nvolseries,replicates,fractionOfVarianceTokeep);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end
disp(['done calculating kmeans.'])


% run makeClusters.m and processVolumeImagingKmeans.m for all folders


rawKmeansOutput=['rawKmeans_' num2str(nkmeans) 'clusters_' num2str(fractionOfVarianceTokeep) 'fractionOfVarianceKept_' num2str(replicates) 'replicates.mat'];
clusterVolFile='clusterVols.mat';
volumeAcquisitionTime=1.2;
nChannels=2;

homeDir=pwd;
for days=1:length(daysToProcess)
    
    cd(daysToProcess{days})
    startDir=pwd;
    display(['processing folder ' startDir])
    currFolders = dir(startDir);
    currFolders=currFolders(3:end);
    clear actuallyAFolder
    for i=1:length(currFolders)
        actuallyAFolder(i)=currFolders(i).isdir;
    end
    currFolders(~actuallyAFolder)=[];
    
    for i=1:length(currFolders)
        cd(currFolders(i).name)
        
        if exist('leftLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['leftLobe\' folderToRunKmeansOn{j}])
                    cd(['leftLobe\' folderToRunKmeansOn{j}])
                    try
                        [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
                        [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
                    catch
                    end
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        if exist('rightLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['rightLobe\'  folderToRunKmeansOn{j}])
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                    try
                        [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
                        [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
                    catch
                    end
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        cd(startDir)
    end
    cd(homeDir)
    drawnow
    close all
end
disp(['done creating clusters and responses.'])
%
% copy clusterVolumes .mat files to dropbox
clear all
folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';

daysToProcess{1}='190211_orcoOX5gcamp_behaviorandimaging_copy';

fileToCopy='clusterResponses_';
destinationFolder='C:\Users\mac0456\Dropbox\flyimaging';

homeDir=pwd;
for days=1:length(daysToProcess)
    cd(daysToProcess{days})
    startDir=pwd;
    display(['processing folder ' startDir])
    currFolders = dir(startDir);
    currFolders=currFolders(3:end);
    clear actuallyAFolder
    for i=1:length(currFolders)
        actuallyAFolder(i)=currFolders(i).isdir;
    end
    currFolders(~actuallyAFolder)=[];
    
    for i=1:length(currFolders)
        cd(currFolders(i).name)
        
        if exist('leftLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['leftLobe\' folderToRunKmeansOn{j}])
                    cd(['leftLobe\' folderToRunKmeansOn{j}])
                    if ~exist([destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\leftLobe'])
                        mkdir([destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\leftLobe'])
                    end
                    [status message]=copyfile([startDir '\' currFolders(i).name '\leftLobe\' folderToRunKmeansOn{j} '\' fileToCopy folderToRunKmeansOn{j} '.mat'],[destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\leftLobe\' fileToCopy folderToRunKmeansOn{j} '.mat']);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        if exist('rightLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['rightLobe\'  folderToRunKmeansOn{j}])
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                    if ~exist([destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\rightLobe'])
                        mkdir([destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\rightLobe'])
                    end
                    [status message]=copyfile([startDir '\' currFolders(i).name '\rightLobe\' folderToRunKmeansOn{j} '\' fileToCopy folderToRunKmeansOn{j} '.mat'],[destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\rightLobe\' fileToCopy folderToRunKmeansOn{j} '.mat']);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end