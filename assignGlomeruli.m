clear all
close all
% assign glomeruli
% Matt Churgin, August 2018
% initialization
leftLobe=0;
minResponse=0.1; % remove clusters with maximum response less than this
odorToPhysWeight=100/100; % should be between 0 and 1
showIntermediateFigs=1;
compressDoORdata=0; % if you want to take log of DooR data to compress (DoOR data is ORN, PNs have compressed signal)
normalizeWithinGlom=0; % if you want to z-score odor responses within each glomerulus

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
pubZ=publishedOR.publishedOR.gh146zCentroid;
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
if normalizeWithinGlom
    for i=1:size(pubOR,2)
        if nanstd(pubOR(:,i))>0
            pubOR(:,i)=(pubOR(:,i)-nanmean(pubOR(:,i)))./nanstd(pubOR(:,i));
        else
            pubOR(:,i)=(pubOR(:,i)-nanmean(pubOR(:,i)));
        end
    end
else
    for i=1:size(pubOR,1)
        if nanstd(pubOR(i,:))>0
            pubOR(i,:)=(pubOR(i,:)-nanmedian(pubOR(i,:)))./max(pubOR(i,:));
        else
            pubOR(i,:)=(pubOR(i,:)-nanmedian(pubOR(i,:)));
        end
    end
end

% centroid data is for the right antenna lobe
if leftLobe
    pubX=-(pubX-nanmedian(pubX))/max(pubX); % flip x axis centroids if looking at a left antenna lobe
else
    pubX=(pubX-nanmedian(pubX))/max(pubX);
end
pubY=(pubY-nanmedian(pubY))/max(pubY);
pubZ=-(pubZ-nanmedian(pubZ))/max(pubZ); % minus sign flips the z centroids to match our data (lower z means more ventral)

% calculate max response over all time points
maxResponse=max(grnResponseNorm,[],3);

myOR=maxResponse(2:13,:); % omit air (odor 1
todelete=max(myOR)<minResponse;
myOR(:,todelete)=NaN;
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
cX(todelete)=NaN;
cY(todelete)=NaN;
cZ(todelete)=NaN;
% cX=(cX-nanmean(cX))/max(cX);
% cY=(cY-nanmean(cY))/max(cY);
% cZ=(cZ-nanmean(cZ))/max(cZ);
cX=(cX-nanmedian(cX))/max(cX);
cY=(cY-nanmedian(cY))/max(cY);
cZ=(cZ-nanmedian(cZ))/max(cZ);


% mean center and normalize each cluster's odor response
% normalize each cluster's odor response (to get relative activation
if normalizeWithinGlom
    for i=1:size(myOR,2)
        if std(myOR(:,i))>0
            myOR(:,i)=(myOR(:,i)-mean(myOR(:,i)))./std(myOR(:,i));
        else
            myOR(:,i)=(myOR(:,i)-mean(myOR(:,i)));
        end
    end
else
    for i=1:size(myOR,1)
        if nanstd(myOR(i,:))>0
            myOR(i,:)=(myOR(i,:)-nanmedian(myOR(i,:)))./max(myOR(i,:));
        else
            myOR(i,:)=(myOR(i,:)-nanmedian(myOR(i,:)));
        end
    end
end

% calculate euclidean distance between each cluster and published responses
% for these odors
odorDist=NaN*zeros(size(myOR,2),size(pubOR,2));
odorDistAcross=NaN*zeros(size(myOR,2),size(pubOR,2));
physDist=NaN*zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
        if any(pubOR(:,j))
            if any(myOR(:,i))
                %odorDist(i,j)=sqrt(nansum((myOR(:,i)-pubOR(:,j)).^2))/sqrt(sum(isfinite(pubOR(:,j))));
                odorDist(i,j)=sqrt(nansum(((myOR(:,i)-pubOR(:,j)))).^2)/sqrt(nansum(((myOR(:,i)+pubOR(:,j)))).^2);
                
                %odorDist(i,j)=sqrt(nansum(((myOR(:,i)-pubOR(:,j))./myORraw(:,i)).^2))./sqrt(sum(isfinite(pubOR(:,j))));
                
                physDist(i,j)=sqrt((cX(i)-pubX(j)).^2+(cY(i)-pubY(j)).^2+(cZ(i)-pubZ(j)).^2);
            end
        end
    end
end

% normalize odorDist and physDist for equal contributions
odorDist=odorDist/max(max(odorDist));
physDist=physDist/max(max(physDist));

compositeDist=(1-odorToPhysWeight)*odorDist+odorToPhysWeight*physDist;

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
    imagesc(odorDist)
    %set(gca,'FontSize',20)
    set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
    xtickangle(30)
    xlabel('Glomeruli','FontSize',20)
    ylabel('Cluster #','FontSize',20)
    title('Odor Panel Response Euclidean Distance','FontSize',20)
    
    figure;
    imagesc(physDist)
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
        text(clusterInfoU{j}.Centroid(1),clusterInfoU{j}.Centroid(2),clusterInfoU{j}.Centroid(3),pubGlomNames{glomMinimizing(j)})
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
    text(clusterInfoU{uniqueClusters(j)}.Centroid(1),clusterInfoU{uniqueClusters(j)}.Centroid(2),clusterInfoU{uniqueClusters(j)}.Centroid(3),pubGlomNames{uniqueGloms(j)})
    hold on
end








