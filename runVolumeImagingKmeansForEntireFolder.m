%% run VolumeImagingkKmeans.m for all flies in a directory
% Matt Churgin, September 2018
clear all

folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';

daysToProcess{1}='181012_pairedbehaviorimaging_gh146';


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

%%
% run makeClusters.m and processVolumeImagingKmeans.m for all folders

folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';
% 
% 
daysToProcess{1}='180831_pairedbehaviorandimaging';
daysToProcess{2}='180906';
daysToProcess{3}='180911_pairedbehaviorimaging';
daysToProcess{4}='180914_pairedbehaviorimaging';
daysToProcess{5}='180918_pairedbehaviorimaging';
daysToProcess{6}='180920_pairebehaviorimaging';
daysToProcess{7}='180925_pairedbehaviorimaging';
daysToProcess{8}='181002_pairedbehaviorimaging_gh146';
daysToProcess{9}='181003_pairedbehaviorimaging_gh146';
daysToProcess{10}='181010_pairedbehaviorimaging_gh146';
daysToProcess{11}='181012_pairedbehaviorimaging_gh146';


rawKmeansOutput='rawKmeans_80clusters_0.75fractionOfVarianceKept_50replicates.mat';
clusterVolFile='clusterVols.mat';
volumeAcquisitionTime=1.2;
nChannels=2;

homeDir=pwd;
for days=1:length(daysToProcess)
    close all
    drawnow
    
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
                    [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
                    [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        if exist('rightLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['rightLobe\'  folderToRunKmeansOn{j}])
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                    [clusterVolU clusterInfoU]=makeClusters(rawKmeansOutput);
                    [grnResponse t] = processVolumeImagingKmeans(clusterVolFile,rawKmeansOutput,volumeAcquisitionTime,nChannels);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end
disp(['done creating clusters and responses.'])

%% copy clusterVolumes .mat files to dropbox
clear all
folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';


daysToProcess{1}='181010_pairedbehaviorimaging_gh146';


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
                        [status message]=copyfile([startDir '\' currFolders(i).name '\leftLobe\' folderToRunKmeansOn{j} '\' fileToCopy folderToRunKmeansOn{j} '.mat'],[destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\leftLobe\' fileToCopy folderToRunKmeansOn{j} '.mat']);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        if exist('rightLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['rightLobe\'  folderToRunKmeansOn{j}])
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                        [status message]=copyfile([startDir '\' currFolders(i).name '\rightLobe\' folderToRunKmeansOn{j} '\' fileToCopy folderToRunKmeansOn{j} '.mat'],[destinationFolder '\' daysToProcess{days} '\' currFolders(i).name '\rightLobe\' fileToCopy folderToRunKmeansOn{j} '.mat']);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end