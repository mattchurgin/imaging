% assign glomeruli
% calculate two odor distances.
% Matt Churgin, August 2018
%% load data and compute spatial and odor response priors

%clear all
close all

% initialization
leftLobe=0;
minResponse=0.1; % remove clusters with maximum response less than this
showIntermediateFigs=0;
compressDoORdata=1; % if you want to take log of DooR data to compress (DoOR data is ORN, PNs have compressed signal)
omitAtlasZSlices=1;

% weights for each prior
% all 1s for equal weighting
physDistWeight=1;
odorWeight=1;
shapeWeight=1;

filename=uigetfile(); % load processed k means .mat file

load(filename)

% load published odor response matrix and centroids
publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
publishedOR=load(publishedOdorPath);
pubOR=publishedOR.publishedOR.gh146response';
pubNames=publishedOR.publishedOR.gh146receptorNames;
pubGlomNames=publishedOR.publishedOR.gh146glomerulusNames;
pubX=publishedOR.publishedOR.gh146xCentroid;
pubY=publishedOR.publishedOR.gh146yCentroid;
pubZ=-publishedOR.publishedOR.gh146zCentroid;% minus sign flips the z centroids to match our data (lower z means more ventral)

% remove z slices
if omitAtlasZSlices
    slicesToRemove=find(pubZ==-5);
    pubZ(slicesToRemove)=NaN;
    pubX(slicesToRemove)=NaN;
    pubY(slicesToRemove)=NaN;
else
    slicesToRemove=[];
end

for j=1:length(pubNames)
    pubNames{j}=num2str(pubNames{j});
    pubGlomNames{j}=num2str(pubGlomNames{j});
end

if compressDoORdata
    % DoOR data is taken from ORNs.  PNs are known to amplify weak signal.
    % Therefore, we compress pubOR signal by taking log of published odor response
    pubOR=log10(pubOR);
    pubOR(isinf(pubOR))=NaN;
    pubOR=pubOR-min(min(pubOR));
end

% mean center and normalize published data

for i=1:size(pubOR,2)
    if nanstd(pubOR(:,i))>0
        pubORWithinGlom(:,i)=(pubOR(:,i)-nanmin(pubOR(:,i)))./(nanmax(pubOR(:,i))-nanmin(pubOR(:,i)));
    else
        pubORWithinGlom(:,i)=(pubOR(:,i)-nanmin(pubOR(:,i)));
    end
end

for i=1:size(pubOR,1)
    if nanstd(pubOR(i,:))>0
        pubORAcrossGlom(i,:)=(pubOR(i,:)-nanmin(pubOR(i,:)))./(max(pubOR(i,:))-min(pubOR(i,:)));
    else
        pubORAcrossGlom(i,:)=(pubOR(i,:)-nanmin(pubOR(i,:)));
    end
end

% calculate max response over all time points
maxResponse=max(grnResponseNorm,[],3);

myOR=maxResponse(2:13,:); % omit air (odor 1
rawClustersToDelete=max(myOR)<minResponse;
myOR(:,rawClustersToDelete)=NaN;
myORraw=myOR;



% centroid data is for the right antenna lobe
if leftLobe
    pubX=-(pubX-nanmin(pubX))/(max(pubX)-min(pubX)); % flip x axis centroids if looking at a left antenna lobe
else
    pubX=(pubX-nanmin(pubX))/(max(pubX)-min(pubX));
end
pubY=(pubY-nanmin(pubY))/(max(pubY)-min(pubY));
pubZ=(pubZ-nanmin(pubZ))/(max(pubZ)-min(pubZ));




% load cluster centroid data and normalize
cX=zeros(1,length(clusterInfoU));
cY=zeros(1,length(clusterInfoU));
cZ=zeros(1,length(clusterInfoU));

for j=1:length(clusterInfoU)
    cX(j)=clusterInfoU{j}.Centroid(1);
    cY(j)=clusterInfoU{j}.Centroid(2);
    cZ(j)=clusterInfoU{j}.Centroid(3);
end
cX(todelete)=NaN;
cY(todelete)=NaN;
cZ(todelete)=NaN;
cX=(cX-nanmin(cX))/(max(cX)-min(cX));
cY=(cY-nanmin(cY))/(max(cY)-min(cY));
cZ=(cZ-nanmin(cZ))/(max(cZ)-min(cZ));

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

% mean center and normalize each cluster's odor response
% normalize each cluster's odor response (to get relative activation

for i=1:size(myOR,2)
    if std(myOR(:,i))>0
        myORWithinGlom(:,i)=(myOR(:,i)-nanmin(myOR(:,i)))./(nanmax(myOR(:,i))-nanmin(myOR(:,i)));
    else
        myORWithinGlom(:,i)=(myOR(:,i)-nanmin(myOR(:,i)));
    end
end

for i=1:size(myOR,1)
    if nanstd(myOR(i,:))>0
        myORAcrossGlom(i,:)=(myOR(i,:)-nanmin(myOR(i,:)))./(nanmax(myOR(i,:))-nanmin(myOR(i,:)));
    else
        myORAcrossGlom(i,:)=(myOR(i,:)-nanmin(myOR(i,:)));
    end
end

% calculate euclidean distance between each cluster and published responses
% for these odors
odorDist=NaN*zeros(size(myOR,2),size(pubOR,2));
odorDistAcross=NaN*zeros(size(myOR,2),size(pubOR,2));
physDist=NaN*zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
                odorDist(i,j)=sqrt(nansum((myORWithinGlom(:,i)-pubORWithinGlom(:,j)).^2))/sqrt(sum(isfinite(pubORWithinGlom(:,j))));
                %odorDist(i,j)=sqrt(nansum(((myOR(:,i)-pubOR(:,j)))).^2)/sqrt(nansum(((myOR(:,i)+pubOR(:,j)))).^2);
                
                %odorDist(i,j)=sqrt(nansum(((myOR(:,i)-pubOR(:,j))./myORraw(:,i)).^2))./sqrt(sum(isfinite(pubOR(:,j))));
                
                
                %odorDistAcross(i,j)=????
                
                physDist(i,j)=sqrt((cX(i)-pubX(j)).^2+(cY(i)-pubY(j)).^2+(cZ(i)-pubZ(j)).^2);
                
                % omit z information
                %physDist(i,j)=sqrt((cX(i)-pubX(j)).^2+(cY(i)-pubY(j)).^2);
    end
end

%% compute glomerulus shape prior
% use 2d cross-correlation between each cluster maximum 2d projection and
% the glomerulus cross-section

% create 2d maximum projections for quantifying shape prior
clusterProj=cell(1,size(myOR,2));
for i=1:size(myOR,2)
    clusterTemp=double(any(clusterVolU{i},3));
    
    % another option: take only the centroid slice to compare to atlas
    %clusterTemp=double(clusterVolU{i}(:,:,round(clusterInfoU{i}.Centroid(3))));
    
    % set pixels outside of cluster equal to -1 to penalize partial
    % overlaps
    clusterTemp(clusterTemp==0)=-1;
    
    clusterProj{i}=clusterTemp;
end

glomProj=cell(1,size(pubOR,2));
for i=1:size(pubOR,2)
    glomTemp=double(publishedOR.publishedOR.gh146glomCentroidMask{i});
    glomTemp(glomTemp==0)=-1;
    
    glomProj{i}=glomTemp;
end

% use maximum projection of ech cluster and convolve with door glomeruluar
% atlas to find cluster that maximizes shape correlation
shapePriorNorm=zeros(size(myOR,2),size(pubOR,2));

xpeaks=zeros(size(myOR,2),size(pubOR,2));
ypeaks=zeros(size(myOR,2),size(pubOR,2));
tic
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
        tempcorrnorm=normxcorr2(clusterProj{i},glomProj{j});
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
        disp(['calculated cross-correlation for cluster ' num2str(i) ', glomerulus ' num2str(j)])
    end
end
disp(['time elapsed to compute cross-correlations: ' num2str(toc) ' seconds'])

%% combine priors and assign glomeruli
% normalize distance matrices between 0 and 1
%odorDistNormed=(odorDist-min(odorDist(:)))/max(odorDist(:)-min(odorDist(:)));
odorDistNormed=odorDist;
physDistNormed=(physDist-min(physDist(:)))/max(physDist(:)-min(physDist(:)));
shapePriorNormed=(shapePriorNorm-min(shapePriorNorm(:)))/max(shapePriorNorm(:)-min(shapePriorNorm(:)));

toFillP=find(any(physDistNormed)==0);
physDistNormed(:,toFillP)=nanmean(physDistNormed(:));
physDistNormed(rawClustersToDelete,:)=NaN;

toFillS=find(any(shapePriorNormed)==0);
shapePriorNormed(:,toFillS)=nanmean(shapePriorNormed(:));
shapePriorNormed(rawClustersToDelete,:)=NaN;

toFillO=find(any(odorDistNormed)==0);
odorDistNormed(:,toFillO)=nanmean(odorDistNormed(:));
odorDistNormed(rawClustersToDelete,:)=NaN;

% rawDist(:,:,1)=aweight*log10(odorDistNormed);
% rawDist(:,:,2)=bweight*log10(physDistNormed);
% rawDist(:,:,3)=cweight*log10(1-shapePriorNormed);
% compositeDist=nansum(rawDist,3);


compositeDist = physDistWeight*log10(physDistNormed) + odorWeight*log10(odorDistNormed) + shapeWeight*log10(1-shapePriorNormed);
compositeDist(:,slicesToRemove)=NaN;

if showIntermediateFigs
    figure
    imagesc(myORraw)
    xlabel('Cluster #')
    ylabel('Odor #')
    title('Cluster Maximum dF/F (Raw)')
    set(gca,'FontSize',20)
    
    figure
    imagesc(myOR)
    xlabel('Cluster #')
    ylabel('Odor #')
    title('Cluster Maximum dF/F (Normalized)')
    set(gca,'FontSize',20)
    
    figure
    imagesc(pubOR)
    set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Odor #','FontSize',20)
    title('Glomerular Maximum Response','FontSize',20)
    
    
    figure;
    imagesc(odorDistNormed)
    %set(gca,'FontSize',20)
    set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Cluster #','FontSize',20)
    title('Odor Panel Response Euclidean Distance','FontSize',20)
    
    figure;
    imagesc(physDistNormed)
    %set(gca,'FontSize',20)
    set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Cluster #','FontSize',20)
    title('Physical Euclidean Distance','FontSize',20)
    
    figure;
    imagesc(shapePriorNormed)
    %set(gca,'FontSize',20)
    set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Cluster #','FontSize',20)
    title('Physical Euclidean Distance','FontSize',20)
    
    figure;
    imagesc((compositeDist))
    %set(gca,'FontSize',20)
    set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Cluster #','FontSize',20)
    title('Composite Distance','FontSize',20)
    
end

% Simple Greedy algorithm with one pass through
% find glomerulus that minimizes multiplied distance to each cluster
glomMinimizing=zeros(1,size(myOR,2));
glomMinimizingMatrix=zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    if any(compositeDist(i,:))
        [val ind]=nanmin(compositeDist(i,:));
        glomMinimizing(i)=ind;
        
        [asdf asdf2]=sort(compositeDist(i,:));
        glomMinimizingMatrix(i,:)=asdf2;
    else
        glomMinimizing(i)=0;
        glomMinimizingMatrix(i,:)=NaN;
    end
end


figure
for j=1:length(glomMinimizing)
    if glomMinimizing(j)>0
        p2=patch(isosurface(clusterVolU{j}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
        isonormals(clusterVolU{j},p2)
        text(clusterInfoU{j}.Centroid(1),clusterInfoU{j}.Centroid(2),clusterInfoU{j}.Centroid(3),pubGlomNames{glomMinimizing(j)},'FontSize',15,'FontWeight','Bold')
        hold on
    end
end

uniqueGloms=unique(glomMinimizing);

% find cluster that minimizes multiplied distance to each glomerulus
uniqueClusters=zeros(1,(length(uniqueGloms)));
for i=2:(length(uniqueGloms))
    currGlom=uniqueGloms(i);
    currClusters=find(glomMinimizing==currGlom);
    [mymin myind]=min(compositeDist(currClusters,currGlom));
    uniqueClusters(i)=currClusters(myind);
end

figure
for j=2:length(uniqueGloms)
    p2=patch(isosurface(clusterVolU{uniqueClusters(j)}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(clusterVolU{uniqueClusters(j)},p2)
    text(clusterInfoU{uniqueClusters(j)}.Centroid(1),clusterInfoU{uniqueClusters(j)}.Centroid(2),clusterInfoU{uniqueClusters(j)}.Centroid(3),pubGlomNames{uniqueGloms(j)},'FontSize',15,'FontWeight','Bold')
    hold on
end

%% More complex greedy algorithm makes multiple passes through
compositeDistTemp=compositeDist;

assignmentThreshold=prctile(prctile(compositeDist,25),25);
%assignmentThreshold=0;

assignmentScore=[];
glomerulusAssignment=[];
clusterAssignment=[];
iters=0;
nassignments=1;
while sum(any(compositeDistTemp))>0
%while iters<3
    % find glomerulus that minimizes multiplied distance to each cluster
    glomMinimizing=zeros(1,size(myOR,2));
    glomMinimizingMatrix=zeros(size(myOR,2),size(pubOR,2));
    for i=1:size(myOR,2)
        if any(compositeDistTemp(i,:))
            [val ind]=nanmin(compositeDistTemp(i,:));
            glomMinimizing(i)=ind;
            
            [asdf asdf2]=sort(compositeDistTemp(i,:));
            glomMinimizingMatrix(i,:)=asdf2;
        else
            glomMinimizing(i)=0;
            glomMinimizingMatrix(i,:)=NaN;
        end
    end
    
    uniqueGloms=unique(glomMinimizing);
    
    % find cluster that minimizes composite distance to each glomerulus
    uniqueClusters=zeros(1,(length(uniqueGloms)));
    for i=2:(length(uniqueGloms))
        currGlom=uniqueGloms(i);
        currClusters=find(glomMinimizing==currGlom);
        [mymin myind]=min(compositeDistTemp(currClusters,currGlom));
        assignmentScore(nassignments)=mymin;
        
        uniqueClusters(i)=currClusters(myind);
        compositeDistTemp(currClusters(myind),:)=NaN;
        compositeDistTemp(:,currGlom)=NaN;
        nassignments=nassignments+1;
    end
    
    glomerulusAssignment=[glomerulusAssignment uniqueGloms];
    clusterAssignment=[clusterAssignment uniqueClusters];
    iters=iters+1;
    disp(['completed ' num2str(iters) ' iterations'])
end

todelete=find(glomerulusAssignment==0);
glomerulusAssignment(todelete)=[];
clusterAssignment(todelete)=[];

% remove high scores (poor fit)
highScores=find(assignmentScore>assignmentThreshold);
glomerulusAssignment(highScores)=[];
clusterAssignment(highScores)=[];

figure
for j=1:length(glomerulusAssignment)
    p2=patch(isosurface(clusterVolU{clusterAssignment(j)}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(clusterVolU{clusterAssignment(j)},p2)
    text(clusterInfoU{clusterAssignment(j)}.Centroid(1),clusterInfoU{clusterAssignment(j)}.Centroid(2),clusterInfoU{clusterAssignment(j)}.Centroid(3),pubGlomNames{glomerulusAssignment(j)},'FontSize',15,'FontWeight','Bold')
    hold on
end