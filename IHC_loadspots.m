clear all 
close all
load(['leftLobe_processed/segmented.mat']) % redo this with full resolution image


folderName='leftLobe_normalresolution_spots';
fileNames=dir(folderName);
fileNames=fileNames(3:end);

temp=imread([folderName '/' fileNames(1).name]);
myspots=zeros(size(temp,1),size(temp,2),length(fileNames));
for j=1:length(fileNames)
    temp=imread([folderName '/' fileNames(j).name]);
    myspots(:,:,j)=temp;
    
    if mod(j,10)==0
        disp(['loaded image ' num2str(j) ' of ' num2str(length(fileNames))])
    end
end
myspots=myspots(:,:,end:-1:1);

myspotsBinary=myspots>0; 

%% get glom spots for downsampled matrix
myspotsBinaryDec=myspotsBinary(1:2:end,1:2:end,1:2:end);

glomSpots=cell(1,length(clusterVols2));
for i=1:length(clusterVols2)
    glomSpots{i}=clusterVols2{i}.*myspotsBinaryDec;
    if mod(i,10)==0
       disp(num2str(i)) 
    end
end

% find number of spots in each glom
glomRegionProps=cell(1,length(clusterVols2));

glomThresh=30;
glomWithSpots=zeros(1,length(clusterVols2));
for i=1:length(clusterVols2)
    lblImgtemp = bwlabeln(glomSpots{i});
    glomRegionProps{i}=regionprops(lblImgtemp);
    
    if length(glomRegionProps{i})>glomThresh
        glomWithSpots(i)=1;
    end
    if mod(i,10)==0
        disp(num2str(i))
    end
end
%% plot segmented gloms and glom spots
disp(['plotting results. found ' num2str(length(clusterVols2)) ' voxel islands'])

figure
view(3);
axis tight
camlight
lighting gouraud

hold on

rng('default')
mycmap=hsv(length(clusterVols2));
mycmap=mycmap(randperm(length(mycmap)),:);

randv=randperm(length(clusterVols2));
for i=randv
    p2=patch(isosurface(clusterVols2{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
    isonormals(clusterVols2{i},p2)

    drawnow
end
disp('done plotting segmented gloms')

figure
view(3);
axis tight
camlight
lighting gouraud

hold on

randv=randperm(length(clusterVols2));
for i=randv
    if glomWithSpots(i)
    p2=patch(isosurface(glomSpots{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',1);
    isonormals(glomSpots{i},p2)

    drawnow
    end
end
disp('done plotting glom spots')