% run VolumeImagingkKmeans.m for all flies in a directory

folderToRunKmeansOn='Volumes';

startDir=pwd;
currFolders = dir(pwd);
currFolders=currFolders(3:end);
for i=1:length(currFolders)
    actuallyAFolder(i)=currFolders(i).isdir;
end
currFolders(~actuallyAFolder)=[];

tic
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
disp(['done.  time elapsed = ' num2str(toc/3600) ' hours.'])
