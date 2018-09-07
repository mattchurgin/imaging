% assign glomeruli v3
% Matt Churgin, August 2018
%% load data and compute spatial and odor response priors

clear all
close all

% colors!
mcolors(1,:)=[0, 0.4470, 0.7410];
mcolors(2,:)=[0.8500, 0.3250, 0.0980];

% initialization
leftLobe=0;
minResponse=0.1; % remove clusters with maximum response less than this
showIntermediateFigs=0;
compressDoORdata=0; % if you want to take log of DooR data to compress (DoOR data is ORN, PNs have compressed signal)
omitAtlasZSlices=1;
myPercentage=1; % for normalizing centroid locations

% weights for each prior
% all 1s for equal weighting
physDistWeight=1;
odorWeight=1;
shapeWeight=1;

filename=uigetfile(); % load processed k means .mat file

load(filename)

% LOAD AND PREPROCESS PUBLISHED ODOR RESPONSE AND SPATIAL DATA
publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
publishedOR=load(publishedOdorPath);
pubOR=publishedOR.publishedOR.gh146response';
pubNames=publishedOR.publishedOR.gh146receptorNames;
pubGlomNames=publishedOR.publishedOR.gh146glomerulusNames;
pubX=publishedOR.publishedOR.gh146xCentroid;
pubY=publishedOR.publishedOR.gh146yCentroid;
pubZ=-publishedOR.publishedOR.gh146zCentroid;% minus sign flips the z centroids to match our data (lower z means more ventral)
pubXborder=publishedOR.publishedOR.gh146glomBorderX;
pubYborder=publishedOR.publishedOR.gh146glomBorderY;
pubZborder=publishedOR.publishedOR.gh146glomBorderZ;
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
if omitAtlasZSlices
    slicesToRemove=find(pubZ==-5);
    pubZ(slicesToRemove)=NaN;
    pubX(slicesToRemove)=NaN;
    pubY(slicesToRemove)=NaN;
else
    slicesToRemove=[];
end
pubX=(pubX-prctile(pubX,myPercentage))/(prctile(pubX,100-myPercentage)-prctile(pubX,myPercentage));
pubY=(pubY-prctile(pubY,myPercentage))/(prctile(pubY,100-myPercentage)-prctile(pubY,myPercentage));
pubZ=(pubZ-prctile(pubZ,myPercentage))/(prctile(pubZ,100-myPercentage)-prctile(pubZ,myPercentage));

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

% normalized within glomerulus response across odors (preserves relative
% odor activation within each glomerulus)
for i=1:size(pubOR,2)
    if nanstd(pubOR(:,i))>0
        pubORWithinGlom(:,i)=(pubOR(:,i)-nanmin(pubOR(:,i)))./(nanmax(pubOR(:,i))-nanmin(pubOR(:,i)));
    else
        pubORWithinGlom(:,i)=(pubOR(:,i)-nanmin(pubOR(:,i)));
    end
end

% normalized across glomerulus response for each odor (preserves relative
% activation across glomeruli)
for i=1:size(pubOR,1)
    if nanstd(pubOR(i,:))>0
        pubORAcrossGlom(i,:)=(pubOR(i,:)-nanmin(pubOR(i,:)))./(max(pubOR(i,:))-min(pubOR(i,:)));
    else
        pubORAcrossGlom(i,:)=(pubOR(i,:)-nanmin(pubOR(i,:)));
    end
end

% LOAD AND PREPROCESS CLUSTER DATA
% calculate max response over all time points
maxResponse=max(grnResponse,[],3);
myOR=maxResponse(2:13,:); % omit air (odor 1
rawClustersToDelete=max(myOR)<minResponse;
myOR(:,rawClustersToDelete)=NaN;
myORraw=myOR;

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

% mean center and normalize each cluster's odor response
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

% create sorted versions of odor response matrices
% can be used to measure rank comparison between classified and published
% odor responses
myORRank=NaN*zeros(size(myORAcrossGlom,1),size(myORAcrossGlom,2));
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

pubORRank=NaN*zeros(size(pubORAcrossGlom,1),size(pubORAcrossGlom,2));
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
end

% calculate spearman (rank) correlation between each cluster and each
% glomerulus
odorRankCorr=NaN*zeros(size(myORRank,2),size(pubORRank,2));
odorRankCorrP=NaN*zeros(size(myORRank,2),size(pubORRank,2));
for i=1:size(myORRank,2)
    for j=1:size(pubORRank,2)
        myORtemp=myORRank(:,i);
        pubORtemp=pubORRank(:,j);
        
        availableOdors=find(isfinite(pubORtemp));
        if any(availableOdors)
            myORavailable=myORtemp(availableOdors);
            pubORavailable=pubORtemp(availableOdors);
            
            % re-rank myORRank to include only available odors
            [vals inds]=sort(myORavailable);
            [myC myP]=corrcoef(pubORavailable,inds);
            odorRankCorr(i,j)=myC(1,2);
            odorRankCorrP(i,j)=myP(1,2);
        end
    end
end

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
% FOR PHYSICAL CENTROID DISTANCE AND ODOR RESPONSE SPACE
odorDist=NaN*zeros(size(myOR,2),size(pubOR,2));
odorDistAcross=NaN*zeros(size(myOR,2),size(pubOR,2));
physDist=NaN*zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
        odorDist(i,j)=sqrt(nansum((myORAcrossGlom(:,i)-pubORAcrossGlom(:,j)).^2))/sqrt(sum(isfinite(pubORAcrossGlom(:,j))));
        physDist(i,j)=sqrt((cX(i)-pubX(j)).^2+(cY(i)-pubY(j)).^2+(cZ(i)-pubZ(j)).^2);
    end
end

%% compute glomerulus shape prior
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
        disp(['calculated cross-correlation for cluster ' num2str(i) ' of ' num2str(size(myOR,2)) ', glomerulus ' num2str(j)])
    end
end
disp(['time elapsed to compute cross-correlations: ' num2str(toc) ' seconds'])

%% combine priors and assign glomeruli
% normalize distance matrices between 0 and 1
odorDistNormed=(odorDist-min(odorDist(:)))/(max(odorDist(:))-min(odorDist(:)));
odorCorrNormed=1-(odorRankCorr-min(odorRankCorr(:)))/(max(odorRankCorr(:))-min(odorRankCorr(:)));
physDistNormed=(physDist-prctile(physDist(:),1))/(max(physDist(:))-min(physDist(:)));

shapePriorN=1-shapePriorNorm; % flip shape prior so lower scores are better
shapePriorNormed=(shapePriorN-prctile(shapePriorN(:),1))/(max(shapePriorN(:))-min(shapePriorN(:)));

% make sure minimum is greater than zero
physDistNormed=physDistNormed-2*min(physDistNormed(:));
shapePriorNormed=shapePriorNormed-2*min(shapePriorNormed(:));

toFillP=find(any(physDistNormed)==0);
physDistNormed(:,toFillP)=nanmean(physDistNormed(:));
physDistNormed(rawClustersToDelete,:)=NaN;

toFillS=find(any(shapePriorNormed)==0);
shapePriorNormed(:,toFillS)=nanmean(shapePriorNormed(:));
shapePriorNormed(rawClustersToDelete,:)=NaN;

toFillO=find(any(odorCorrNormed)==0);
odorCorrNormed(:,toFillO)=nanmean(odorCorrNormed(:));
odorCorrNormed(rawClustersToDelete,:)=NaN;

odorCorrNormed(:,slicesToRemove)=NaN;
physDistNormed(:,slicesToRemove)=NaN;
shapePriorNormed(:,slicesToRemove)=NaN;


% sort distance matrices along rows (find glom rank for each cluster)
[garbage physClusterRankTemp]=sort(physDistNormed,2);
for i=1:size(physDistNormed,1)
    for j=1:size(physDistNormed,2)
        physClusterRank(i,j)=find(physClusterRankTemp(i,:)==j);
    end
end
[garbage shapeClusterRankTemp]=sort(shapePriorNormed,2);
for i=1:size(physDistNormed,1)
    for j=1:size(physDistNormed,2)
        shapeClusterRank(i,j)=find(shapeClusterRankTemp(i,:)==j);
    end
end
[garbage odorClusterRankTemp]=sort(odorCorrNormed,2);
for i=1:size(physDistNormed,1)
    for j=1:size(physDistNormed,2)
        odorClusterRank(i,j)=find(odorClusterRankTemp(i,:)==j);
    end
end


odorClusterRank(:,slicesToRemove)=NaN;
shapeClusterRank(:,slicesToRemove)=NaN;
physClusterRank(:,slicesToRemove)=NaN;

odorClusterRank(rawClustersToDelete,:)=NaN;
shapeClusterRank(rawClustersToDelete,:)=NaN;
physClusterRank(rawClustersToDelete,:)=NaN;

if showIntermediateFigs
    figure
    imagesc(myOR)
    xlabel('Cluster #')
    ylabel('Odor #')
    title('Cluster Maximum dF/F (Normalized)')
    set(gca,'FontSize',20)
    
    figure
    imagesc(pubOR)
    set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Odor #','FontSize',20)
    title('Glomerular Maximum Response','FontSize',20)
end

figure;
subplot(3,1,1)
imagesc(odorCorrNormed)
set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Odor Panel Rank Correlation (Normed)','FontSize',20)

subplot(3,1,2)
imagesc(physDistNormed)
set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Physical Euclidean Distance (Normed)','FontSize',20)

subplot(3,1,3)
imagesc(shapePriorNormed)
set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Shape Score (Normed)','FontSize',20)


figure;
subplot(3,1,1)
imagesc(odorClusterRank)
set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Odor Glom Rank','FontSize',20)

subplot(3,1,2)
imagesc(physClusterRank)
set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Physical Euclidean Distance Glom Rank','FontSize',20)

subplot(3,1,3)
imagesc(shapeClusterRank)
set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Shape Glom Rank','FontSize',20)

figure
for j=1:length(clusterVolU)
    p2=patch(isosurface(clusterVolU{j}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(clusterVolU{j},p2)
    text(clusterInfoU{j}.Centroid(1),clusterInfoU{j}.Centroid(2),clusterInfoU{j}.Centroid(3),[num2str(j)],'FontSize',15,'FontWeight','Bold')
    hold on
end
%% Use random permutations to try and improve greedy algorithm
nshuffles=1;
% what if you use cluster ranks as distance input?  also, do the thing you
% said in your email, i.e., create three composite dists and try to
% optimize each, then choose the one with the lowest total score based on
% rank.......
[output] = assignGloms(nshuffles,odorCorrNormed,physDistNormed,shapePriorNormed,slicesToRemove,odorRankCorr);

figure
plot(output.optimizedScore)
hold on
plot(1:length(output.optimizedScore),output.optimizedScore(1)*ones(1,length(output.optimizedScore)),'k--','LineWidth',2)
plot(output.bestIndex,output.optimizedScore(output.bestIndex),'ro','LineWidth',2)
xlabel('Trial #')
ylabel('Total Assignment Score')
legend('Trials','Vanilla Greedy','Winner')
legend boxoff
box off
set(gca,'FontSize',15)

figure;
plot(output.optimizedScore,output.totalOdorScore,'o')
hold on
plot(output.optimizedScore,output.totalOdorScoreShuffled,'r.')
plot(output.optimizedScore(1),output.totalOdorScoreShuffled(1),'ko','LineWidth',2)
legend('Unshuffled','Shuffled')
xlabel('Total Assignment Score')
ylabel('Total Odor Score')
legend boxoff
box off
set(gca,'FontSize',15)

figure
for j=1:length(output.glomerulusAssignmentFinal)
    p2=patch(isosurface(clusterVolU{output.clusterAssignmentFinal(j)}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(clusterVolU{output.clusterAssignmentFinal(j)},p2)
    text(clusterInfoU{output.clusterAssignmentFinal(j)}.Centroid(1),clusterInfoU{output.clusterAssignmentFinal(j)}.Centroid(2),clusterInfoU{output.clusterAssignmentFinal(j)}.Centroid(3),pubGlomNames{output.glomerulusAssignmentFinal(j)},'FontSize',15,'FontWeight','Bold')
    hold on
end

%save('assignedGlomeruli.mat','output','compositeDist','myORAcrossGlom','pubORAcrossGlom','odorDistNormed','physDistNormed','shapePriorNormed','pubGlomNames')