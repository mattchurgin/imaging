% assign glomeruli v3
% Matt Churgin, August 2018
% load data and compute spatial and odor response distance and rank correlation matrices

clear all
close all

% colors for plotting
mcolors(1,:)=[0, 0.4470, 0.7410];
mcolors(2,:)=[0.8500, 0.3250, 0.0980];

% initialization
%orn=input('ORN (0) or PN (1)? ');
orn=1;

% assuming data is in a folder named leftLobe or rightLobe
% automatically extract which lobe current data is from
currfolder=pwd;
splits=strfind(currfolder,'/');
lobename=currfolder((splits(end)+1):end);
leftLobe=strfind(lobename,'left');
leftLobe=1-length(leftLobe); % do this because in actuality lobes labelled as left are in fact on the right
minResponse=0.05; % remove clusters with maximum response less than this
myPercentage=1; % for normalizing centroid locations
minOdors=4; % omit odor rank correlations when less than min odors are available in DOOR dataset

clusterfilename=uigetfile(); % load file with clusters and responses
savesuffix=clusterfilename(18:(end-4));

load(clusterfilename)
showClusters(clusterVolU,clusterInfoU)

clustersManuallyOmitted=input('Any clusters to omit? ');
omitAtlasZslices=input('omit any atlas z slices?  (5 is ventral, 1 is dorsal) : ');

% LOAD AND PREPROCESS PUBLISHED ODOR RESPONSE AND SPATIAL DATA
publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
publishedOR=load(publishedOdorPath);
if orn==0
    
elseif orn==1
    pubOR=publishedOR.publishedOR.gh146response';
    pubNames=publishedOR.publishedOR.gh146receptorNames;
    pubGlomNames=publishedOR.publishedOR.gh146glomerulusNames;
    pubX=publishedOR.publishedOR.gh146xCentroid;
    pubY=publishedOR.publishedOR.gh146yCentroid;
    pubZ=-publishedOR.publishedOR.gh146zCentroid;% minus sign flips the z centroids to match our data (lower z means more ventral)
    pubXborder=publishedOR.publishedOR.gh146glomBorderX;
    pubYborder=publishedOR.publishedOR.gh146glomBorderY;
    pubZborder=publishedOR.publishedOR.gh146glomBorderZ;
end
pubXborderAll=[];
pubYborderAll=[];
for j=1:length(pubXborder)
    pubXborderAll=[pubXborderAll pubXborder{j}];
    pubYborderAll=[pubYborderAll pubYborder{j}];
end
xscale=size(clusterVolU{1},2)/max(pubXborderAll); % find x and y scaling from atlas to 2-photon images
yscale=size(clusterVolU{1},1)/max(pubYborderAll);

% flip x axis if data is a left lobe
if leftLobe
    pubX=-pubX;
end

% remove z slices
slicesToRemove=[];
for i=1:length(omitAtlasZslices)
    slicesToRemove=[slicesToRemove find(pubZ==-omitAtlasZslices(i))];
    pubZ(slicesToRemove)=NaN;
    pubX(slicesToRemove)=NaN;
    pubY(slicesToRemove)=NaN;
end
pubX=(pubX-prctile(pubX,myPercentage))/(prctile(pubX,100-myPercentage)-prctile(pubX,myPercentage));
pubY=(pubY-prctile(pubY,myPercentage))/(prctile(pubY,100-myPercentage)-prctile(pubY,myPercentage));
pubZ=(pubZ-prctile(pubZ,myPercentage))/(prctile(pubZ,100-myPercentage)-prctile(pubZ,myPercentage));

% find pair-wise distances between all published glomeruli locations
pubPairwiseDist=NaN*zeros(length(pubX),length(pubX));
for i=1:length(pubX)
    for j=1:length(pubX)
        if i~=j
        pubPairwiseDist(i,j)=sqrt((pubX(i)-pubX(j)).^2+(pubY(i)-pubY(j)).^2+(pubZ(i)-pubZ(j)).^2);
        end
    end
end

% convert pairwise distances into rank distances
pubPairwiseDistRank=NaN*zeros(length(pubX),length(pubX));
for i=1:length(pubX)
    [vals inds]=sort(pubPairwiseDist(i,:));
    % remove nans
    todelete=find(isnan(vals));
    vals(todelete)=[];
    inds(todelete)=[];
    
    % save rank order
    for j=1:length(inds)
        pubPairwiseDistRank(i,inds(j))=j;
    end
end

% get published glomerulus names
for j=1:length(pubNames)
    pubNames{j}=num2str(pubNames{j});
    pubGlomNames{j}=num2str(pubGlomNames{j});
end

% LOAD AND PREPROCESS CLUSTER DATA
% calculate max response over all time points

% time points in which to calculate average response for summarizing odor ranks for each cluster
% these indices are hard-coded for now (time points 6 s - 12 s)
timepointstoconsider=[5:10]; 

maxResponse=max(grnResponse,[],3);
maxResponse=median(grnResponse(:,:,timepointstoconsider),3); 

% grnResponsePositive=grnResponse;
% grnResponsePositive(grnResponsePositive<0)=0;
% sumResponse=sum(grnResponsePositive,3);

myOR=maxResponse(2:13,:); % omit air (odor 1)
myORraw=myOR;
rawClustersToDelete=max(abs(myOR))<minResponse;
%rawClustersToDelete=logical(zeros(1,size(myOR,2)));
rawClustersToDelete(clustersManuallyOmitted)=1;
myOR(:,rawClustersToDelete)=NaN;


% load cluster centroid data and normalize
cX=zeros(1,length(clusterInfoU));
cY=zeros(1,length(clusterInfoU));
cZ=zeros(1,length(clusterInfoU));

for j=1:length(clusterInfoU)
    cX(j)=clusterInfoU{j}.Centroid(1);
    cY(j)=clusterInfoU{j}.Centroid(2);
    cZ(j)=clusterInfoU{j}.Centroid(3);
end
cX(rawClustersToDelete)=NaN;
cY(rawClustersToDelete)=NaN;
cZ(rawClustersToDelete)=NaN;
cX=(cX-prctile(cX,myPercentage))/(prctile(cX,100-myPercentage)-prctile(cX,myPercentage));
cY=(cY-prctile(cY,myPercentage))/(prctile(cY,100-myPercentage)-prctile(cY,myPercentage));
cZ=(cZ-prctile(cZ,myPercentage))/(prctile(cZ,100-myPercentage)-prctile(cZ,myPercentage));

% create sorted versions of odor response matrices
% can be used to measure rank comparison between classified and published
% odor responses
myORRank=NaN*zeros(size(myOR,1),size(myOR,2));
for i=1:size(myOR,2)
    [vals inds]=sort(myOR(:,i));
    % remove nans
    todelete=find(isnan(vals));
    vals(todelete)=[];
    inds(todelete)=[];
    
    % save rank order
    for j=1:length(inds)
        myORRank(inds(j),i)=j;
    end
end

pubORRank=NaN*zeros(size(pubOR,1),size(pubOR,2));
for i=1:size(pubOR,2)
    [vals inds]=sort(pubOR(:,i));
    % remove nans
    todelete=find(isnan(vals));
    vals(todelete)=[];
    inds(todelete)=[];
    
    % save rank order
    for j=1:length(inds)
        pubORRank(inds(j),i)=j;
    end
    
    numOdorsAvailablePerGlom(i)=sum(isfinite(pubORRank(:,i)));
end

% calculate spearman (rank) correlation between each cluster and each
% glomerulus
odorRankCorr=NaN*zeros(size(myORRank,2),size(pubORRank,2));
for i=1:size(myORRank,2)
    for j=1:size(pubORRank,2)
        myORtemp=myORRank(:,i);
        pubORtemp=pubORRank(:,j);
        
        availableOdors=find(isfinite(pubORtemp));
        
        if any(availableOdors) && any(myORtemp)
            myORavailable=myORtemp(availableOdors);
            pubORavailable=pubORtemp(availableOdors);
            
            % re-rank myORRank to include only available odors
            [vals inds]=sort(myORavailable);
            [myC myP]=corrcoef(pubORavailable,inds);
            odorRankCorr(i,j)=myC(1,2);
        end
    end
end

fullOdorRankCorr=odorRankCorr;
glomsWithLowOdorData=numOdorsAvailablePerGlom<minOdors;
odorRankCorr(:,glomsWithLowOdorData)=NaN;

% visualize centroids
figure
plot3(pubX,pubY,pubZ,'o')
hold on
for j=1:length(pubX)
    try
        text(pubX(j),pubY(j),pubZ(j), pubGlomNames{j},'FontSize',15,'FontWeight','Bold')
    catch
    end
end
plot3(cX,cY,cZ,'ro')

% CALCULATE EUCLIDEAN DISTANCE BETWEEN EACH CLUSTER AND PUBLISHED DATA
% FOR PHYSICAL CENTROID DISTANCE 
%odorDist=NaN*zeros(size(myOR,2),size(pubOR,2));
physDist=NaN*zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
        physDist(i,j)=sqrt((cX(i)-pubX(j)).^2+(cY(i)-pubY(j)).^2+(cZ(i)-pubZ(j)).^2);
    end
end

% compute glomerulus shape prior
% use 2d cross-correlation between each cluster maximum 2d projection and
% the glomerulus cross-section
nonoverlappenalty=1; % value to set mask outside of cluster boundary

% create 2d maximum projections for quantifying shape prior
clusterProj=cell(1,size(myOR,2));
for i=1:size(myOR,2)
    clusterTemp=double(any(clusterVolU{i},3));
    
    % set pixels outside of cluster equal to -1 to penalize partial
    % overlaps
    clusterTemp(clusterTemp==0)=-nonoverlappenalty;
    
    clusterProj{i}=clusterTemp;
end

% make glomProj from atlas glom borders
glomCentroidMask=cell(1,size(pubOR,2));
glomProj=cell(1,size(pubOR,2));
for i=1:size(pubOR,2)
    % scale coordinates according to size of 2-photon image
    pubXborderScaled{i}=xscale*pubXborder{i};
    pubYborderScaled{i}=yscale*pubYborder{i};
    
    %  convert boundary into logical matrices
    [xx yy]=meshgrid(1:size(clusterVolU{1},2),1:size(clusterVolU{1},1));
    glomCentroidMask{i} = inpolygon(xx,yy,pubXborderScaled{i},pubYborderScaled{i});
    
    glomTemp=double(glomCentroidMask{i});
    glomTemp(glomTemp==0)=-nonoverlappenalty;
    if leftLobe
        glomProj{i}=glomTemp(:,end:-1:1);
    else
        glomProj{i}=glomTemp;
    end
end

% use maximum projection of ech cluster and perform cross-correlation with door glomeruluar
% atlas to find cluster that maximizes shape correlation
shapePriorNorm=zeros(size(myOR,2),size(pubOR,2));

xpeaks=zeros(size(myOR,2),size(pubOR,2));
ypeaks=zeros(size(myOR,2),size(pubOR,2));
tic
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
        %tempcorrnorm=normxcorr2(clusterProj{i},glomProj{j});
        tempcorrnorm=[1 2 1]; % don't run actually 2-d cross corr because i'm not using that shape information for now
        shapePriorNorm(i,j)=max(tempcorrnorm(:));
        
        % save location of maximum cross correlation
        [ypeak xpeak]=find(tempcorrnorm==max(tempcorrnorm(:)));
        try
            xpeaks(i,j)=xpeak;
            ypeaks(i,j)=ypeak;
        catch
            xpeaks(i,j)=xpeak(1);
            ypeaks(i,j)=ypeak(1);
        end
        disp(['calculated cross-correlation for cluster ' num2str(i) ' of ' num2str(size(myOR,2)) ', glomerulus ' num2str(j)])
    end
end

disp(['time elapsed to compute cross-correlations: ' num2str(toc) ' seconds'])

save(['intermediateAssignment_' savesuffix '.mat'], 'myOR', 'pubOR', 'shapePriorNorm', 'physDist', 'fullOdorRankCorr', 'odorRankCorr','pubPairwiseDistRank', 'rawClustersToDelete','slicesToRemove', 'pubNames', 'pubGlomNames','clusterfilename','clustersManuallyOmitted','omitAtlasZslices','minResponse','savesuffix','cX','cY','cZ','pubX','pubY','pubZ')

%% run assign glomeruli
clear all
nShuffles=10000;
intermediatefilename='intermediateAssignment_Volumes.mat';
[output besti]=runAssignGloms(intermediatefilename,nShuffles);

%% run assign glomeruli for all folders 

nShuffles=100000;
odorpanelVolume{1}='Volumes';
odorpanelVolume{2}='Volumes2';

daysToProcess{1}='180831_pairedbehaviorandimaging';
daysToProcess{2}='180906';
daysToProcess{3}='180911_pairedbehaviorimaging';
daysToProcess{4}='180914_pairedbehaviorimaging';
daysToProcess{5}='180918_pairedbehaviorimaging';
daysToProcess{6}='180925_pairedbehaviorimaging';
daysToProcess{7}='181002_pairedbehaviorimaging_gh146';
daysToProcess{8}='181003_pairedbehaviorimaging_gh146';

fileToLoad='intermediateAssignment_';

lobes{1}='leftLobe';
lobes{2}='rightLobe';

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
        for currLobe=1:length(lobes)
            if exist(lobes{currLobe})
                cd(lobes{currLobe})
                for j=1:length(odorpanelVolume)
                    intermediatefilename=[fileToLoad odorpanelVolume{j} '.mat'];
                    if exist(intermediatefilename)
                        [output besti]=runAssignGloms(intermediatefilename,nShuffles);
                        drawnow
                        close all
                    end
                end
                cd([startDir '\' currFolders(i).name])
            end
        end
        cd(startDir)
    end
    cd(homeDir)
end
