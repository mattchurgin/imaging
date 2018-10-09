function [output besti] = runAssignGloms(intermediatefilename,nshuffles)

load(intermediatefilename)
load(clusterfilename)

showIntermediateFigs=0;
% normalize distance matrices 
%odorDistNormed=(odorDist-min(odorDist(:)))/(max(odorDist(:))-min(odorDist(:)));
odorCorrNormed=1-(odorRankCorr+1)/2;
physDistNormed=physDist/max(physDist(:));

shapePriorN=1-shapePriorNorm; % flip shape prior so lower scores are better
shapePriorNormed=(shapePriorN)/(max(shapePriorN(:)));

% make sure minimum is greater than zero
%physDistNormed=physDistNormed-2*min(physDistNormed(:));
%shapePriorNormed=shapePriorNormed-2*min(shapePriorNormed(:));

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
% 
% figure;
% subplot(3,1,1)
% imagesc(odorCorrNormed)
% set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
% xtickangle(30)
% xlabel('Glomeruli','FontSize',20)
% ylabel('Cluster #','FontSize',20)
% title('Odor Panel Rank Correlation (Normed)','FontSize',20)
% 
% subplot(3,1,2)
% imagesc(physDistNormed)
% set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
% xtickangle(30)
% xlabel('Glomeruli','FontSize',20)
% ylabel('Cluster #','FontSize',20)
% title('Physical Euclidean Distance (Normed)','FontSize',20)
% 
% subplot(3,1,3)
% imagesc(shapePriorNormed)
% set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
% xtickangle(30)
% xlabel('Glomeruli','FontSize',20)
% ylabel('Cluster #','FontSize',20)
% title('Shape Score (Normed)','FontSize',20)
% 
% 
% figure;
% subplot(3,1,1)
% imagesc(odorClusterRank)
% set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
% xtickangle(30)
% xlabel('Glomeruli','FontSize',20)
% ylabel('Cluster #','FontSize',20)
% title('Odor Glom Rank','FontSize',20)
% 
% subplot(3,1,2)
% imagesc(physClusterRank)
% set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
% xtickangle(30)
% xlabel('Glomeruli','FontSize',20)
% ylabel('Cluster #','FontSize',20)
% title('Physical Euclidean Distance Glom Rank','FontSize',20)
% 
% subplot(3,1,3)
% imagesc(shapeClusterRank)
% set(gca,'xtick',1:length(pubGlomNames),'xticklabel',string(pubGlomNames),'FontSize',10)
% xtickangle(30)
% xlabel('Glomeruli','FontSize',20)
% ylabel('Cluster #','FontSize',20)
% title('Shape Glom Rank','FontSize',20)

% run classification
close all
clusterCentroids=[cX' cY' cZ'];
glomCentroids=[pubX' pubY' pubZ'];

[output] = assignGloms(nshuffles,odorCorrNormed,physDistNormed,shapePriorNormed,slicesToRemove,odorRankCorr,clusterCentroids,glomCentroids);

% calculate final score for classifying and save
distWeight=0.2;
odorWeight=0.0;
pairwiseWeight=0.5;
% calculate score for final classification
finalScoreForClassifying=zeros(1,length(output.totalDistScore));
for i=1:length(output.totalDistScore)
   finalScoreForClassifying(i)=nanmedian(distWeight*output.distScoreHistory{i}-pairwiseWeight*output.pairwiseDistScoreHistory{i}-odorWeight*output.odorScoreHistory{i});
end

[bestv besti]=min(finalScoreForClassifying);
clusterAssignment=output.clusterAssignmentHistory{besti};
glomerulusAssignment=output.glomerulusAssignmentHistory{besti};
%todelete=find(output.optimizedScoreHistory{besti}>output.assignmentThreshold);
%clusterAssignment(todelete)=[];
%glomerulusAssignment(todelete)=[];

figure
plot(output.distScoreHistory{besti},output.odorScoreHistory{besti},'o','LineWidth',2)
hold on
plot(output.distScoreHistory{besti},output.odorScoreShuffledHistory{besti},'rx','LineWidth',2)
legend('Unshuffled','Shuffled')
legend boxoff
xlabel('Normalized Euclidean Distance')
ylabel('Odor Rank Correlation')
box off
set(gca,'FontSize',15)

figure;
plot(output.optimizedScore,output.totalOdorScore,'o')
hold on
plot(output.optimizedScore,output.totalOdorScoreShuffled,'r.')
plot(output.optimizedScore(besti),output.totalOdorScore(besti),'ko','LineWidth',3)
legend('Unshuffled','Shuffled')
xlabel('Total Assignment Score')
ylabel('Median Odor Rank Correlation')
legend boxoff
box off
set(gca,'FontSize',15)

for j=1:length(glomerulusAssignment)
    clusterVolAssigned{j}=clusterVolU{clusterAssignment(j)};
    clusterInfoAssigned{j}=clusterInfoU{clusterAssignment(j)};
    clusterLabels{j}=pubGlomNames{glomerulusAssignment(j)};
end

showClusters(clusterVolAssigned,clusterInfoAssigned,clusterLabels);

% outputs to save
totalOdorScore=output.totalOdorScore(besti);
totalOptimizedScore=output.optimizedScore(besti);
totalDistScore=output.totalDistScore(besti);
totalPairwiseScore=output.totalPairwiseDistScore(besti);
odorScore=output.odorScoreHistory{besti};
distScore=output.distScoreHistory{besti};
pairwiseScore=output.pairwiseDistScoreHistory{besti};
clusterAssignments=output.clusterAssignmentHistory{besti};

save(['assignedClusters_d' num2str(distWeight) '_p' num2str(pairwiseWeight) '_o' num2str(odorWeight) '_' savesuffix '.mat'],'totalOdorScore','totalOptimizedScore','totalDistScore','totalPairwiseScore','odorScore','distScore','pairwiseScore','clusterVolAssigned','clusterInfoAssigned','clusterLabels','besti','odorWeight','distWeight','pairwiseWeight','clusterAssignments')