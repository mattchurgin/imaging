%% run VolumeImagingkKmeans.m for all flies in a directory
% Matt Churgin, September 2018
clear all

folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';

daysToProcess{1}='181004_orcogal4';


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

% daysToProcess{1}='180831_pairedbehaviorandimaging';
% daysToProcess{2}='180906';
% daysToProcess{3}='180911_pairedbehaviorimaging';
% daysToProcess{4}='180914_pairedbehaviorimaging';
% daysToProcess{1}='180918_pairedbehaviorimaging';
% daysToProcess{2}='180925_pairedbehaviorimaging';
% daysToProcess{3}='181002_pairedbehaviorimaging_gh146';
% daysToProcess{4}='181003_pairedbehaviorimaging_gh146';
daysToProcess{1}='180926_orcogal4xuasgcamp';
daysToProcess{2}='181004_orcogal4';
daysToProcess{3}='180925_orcogal4xuasgcamp';

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

folderToRunKmeansOn{1}='Volumes';
folderToRunKmeansOn{2}='Volumes2';

daysToProcess{1}='180831_pairedbehaviorandimaging';
daysToProcess{2}='180906';
daysToProcess{3}='180911_pairedbehaviorimaging';
daysToProcess{4}='180914_pairedbehaviorimaging';
daysToProcess{5}='180918_pairedbehaviorimaging';
daysToProcess{6}='180925_pairedbehaviorimaging';
daysToProcess{7}='181002_pairedbehaviorimaging_gh146';
daysToProcess{8}='181003_pairedbehaviorimaging_gh146';

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
                        status=copyfile([homeDir '\' startDir '\' currFolders(i).name '\' 'leftLobe' fileToCopy '_' folderToRunKmeansOn{j} '.mat'],[destinationFolder '\' currFolders(i).name '\' 'leftLobe' fileToCopy '_' folderToRunKmeansOn{j} '.mat']);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        if exist('rightLobe')
            for j=1:length(folderToRunKmeansOn)
                if exist(['rightLobe\'  folderToRunKmeansOn{j}])
                    cd(['rightLobe\'  folderToRunKmeansOn{j}])
                        status=copyfile([homeDir '\' startDir '\' currFolders(i).name '\' 'rightLobe' fileToCopy '_' folderToRunKmeansOn{j} '.mat'],[destinationFolder '\' currFolders(i).name '\' 'rightLobe' fileToCopy '_' folderToRunKmeansOn{j} '.mat']);
                    cd([startDir '\' currFolders(i).name])
                end
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end