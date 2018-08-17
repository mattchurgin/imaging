
% assign glomeruli
filename=uigetfile(); % load processed k means .mat file
load(filename)

% small change
% small change second computer
% third change
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

% mean center and normalize published data
for i=1:size(pubOR,2)
    if nanstd(pubOR(:,i))>0
        pubOR(:,i)=(pubOR(:,i)-nanmean(pubOR(:,i)))./nanstd(pubOR(:,i));
    else
        pubOR(:,i)=(pubOR(:,i)-nanmean(pubOR(:,i)));
    end
end
pubX=(pubX-nanmean(pubX))/max(pubX);
pubY=(pubY-nanmean(pubY))/max(pubY);
pubZ=(pubZ-nanmean(pubZ))/max(pubZ);

% calculate max response over all time points
maxResponse=max(grnResponseNorm,[],3);

% omit air
myOR=maxResponse(2:13,:);
minResponse=0.1; % remove clusters with maximum response less than this
todelete=max(myOR)<minResponse;
myOR(:,todelete)=NaN;

% mean center and normalize each cluster
for i=1:size(myOR,2)
    if std(myOR(:,i))>0
        myOR(:,i)=(myOR(:,i)-mean(myOR(:,i)))./std(myOR(:,i));
    else
        myOR(:,i)=(myOR(:,i)-mean(myOR(:,i)));
    end
end

% calculate euclidean distance between each cluster and published responses
% for these odors
distM=NaN*zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    for j=1:size(pubOR,2)
        if any(pubOR(:,j))
            if any(myOR(:,i))
                distM(i,j)=sqrt(nansum((myOR(:,i)-pubOR(:,j)).^2))/sqrt(sum(isfinite(pubOR(:,j))));
            end
        end
    end
end

figure
imagesc(myOR)
xlabel('Cluster #')
ylabel('Odor #')
title('Cluster Maximum dF/F')
set(gca,'FontSize',20)

figure
imagesc(pubOR)
set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Odor #','FontSize',20)
title('Glomerular Maximum Response','FontSize',20)


figure;
imagesc(distM)
%set(gca,'FontSize',20)
set(gca,'xtick',1:length(pubNames),'xticklabel',string(pubNames),'FontSize',5)
xtickangle(30)
xlabel('Glomeruli','FontSize',20)
ylabel('Cluster #','FontSize',20)
title('Odor Panel Response Euclidean Distance','FontSize',20)

% find glomerulus that minimizes distance to each cluster
glomMinimizing=zeros(1,size(myOR,2));
glomMinimizingMatrix=zeros(size(myOR,2),size(pubOR,2));
for i=1:size(myOR,2)
    if any(distM(i,:))
        [val ind]=nanmin(distM(i,:));
        glomMinimizing(i)=ind;
        
        [asdf asdf2]=sort(distM(i,:));
        glomMinimizingMatrix(i,:)=asdf2;
    else
        glomMinimizing(i)=0;
        glomMinimizingMatrix(i,:)=NaN;
    end
end

distMtemp=distM;
% find cluster that minimizes distance to each glomerulus
clusterMinimizing=zeros(1,size(pubOR,2));
for i=1:length(glomMinimizing)
    if glomMinimizing>0

        
    end
end

%%
% find cluster that minimizes distance to each glomeruli
for i=1:size(pubOR,2)
    if any(distM(:,i))
        [val ind]=nanmin(distM(:,i));
        clusterMinimizing(i)=ind;
    else
        clusterMinimizing(i)=0;
    end
end

pcl=unique(clusterMinimizing);

for j=2:length(pcl)
    [val ind]=min(distM(pcl(j),:));
    glomMinimizing(pcl(j))=ind;
end
%
% assign glomeruli by finding minimal distance from each cluster
for i=1:size(myOR,2)
    if sum(isfinite(distM(i,:)))>0
        [val ind]=nanmin(distM(i,:));
        clusterName{i}=pubNames{ind};
    else
        clusterName{i}='censored';
    end
end

% next steps:
% need to optimally assign glomeruli according to distances
%% try to assign each cluster to a glomerulus
alreadyTaken=[];
distMcopy=distM;
for i=1:size(pubOR,2)
    if any(distM(:,i))
        confirmi=0;
        ind=[];
        while confirmi~=i
            % find cluster with minimal distance to ith glom
            [val ind]=nanmin(distM(:,i));
            
            % find glom with minimal distance to ind cluster
            [valc confirmi]=nanmin(distM(ind,:));
        end
        
        clusterMinimizing(i)=ind;
        
    else
        clusterMinimizing(i)=0;
    end
end
%%
figure
view(3);
axis tight
camlight
lighting gouraud
hold on
for j=2:length(pcl)
    p2=patch(isosurface(clusterVols{pcl(j)}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.25);
    isonormals(clusterVols{pcl(j)},p2)
end

% for j=1:length(clusterVols)
%     if ~todelete(j)
%         p2=patch(isosurface(clusterVols{j}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.5);
%         isonormals(clusterVols{j},p2)
%     end
%     %drawnow
%     %pause
% end