% compare manual and automated cluster labels
clear all
close all

manualLabelHome='/Users/mattchurgin/Dropbox/flyimaging/analysis/AutomatedGlomerulusAssignmentValidation_GH146/manuallyLabelled';
automatedLabelHome='/Users/mattchurgin/Dropbox/flyimaging';

automatedLabelFilePrefix='assignedClusters_reprocessed_d0.2_p0.5_o0_';

publishedOdorPath='/Users/mattchurgin/Desktop/dblab/mattFunctions/odorpanelResponsesDoOR/odorPanel_12/odorPanel_12_DoORData.mat';
load(publishedOdorPath);

manualLabelledFolders=dir(manualLabelHome);
manualLabelledFolders=manualLabelledFolders(3:end);

automatedLabelFolders=dir(automatedLabelHome);
automatedLabelFolders=automatedLabelFolders(3:end);

labelComparison=cell(1,length(manualLabelledFolders));
individualConfusionMatrix=cell(1,length(manualLabelledFolders));
totalConfusion=zeros(length(publishedOR.gh146glomerulusNames));
for i=1:length(manualLabelledFolders)
    currname=manualLabelledFolders(i).name;
    if strcmp(currname(end),'1')
        currVol='Volumes';
    else
        currVol='Volumes2';
    end
    if strcmp(currname(end-2),'r')
        currLobe='rightLobe';
    else
        currLobe='leftLobe';
    end
    currDate=currname(1:6);
    underscores=strfind(currname,'_');
    currFlyNum=currname((underscores(1)+4):(underscores(2)-1));
    
    currFiles=dir([manualLabelHome '/' manualLabelledFolders(i).name]);
    for j=1:length(currFiles)
       if strfind(currFiles(j).name,'manualLabels')
          manualLabelFile=currFiles(j).name;
          % load manual label file
          load([manualLabelHome '/' manualLabelledFolders(i).name '/' manualLabelFile])
           break
       end
    end
    
    manualClusterLabels=clusterLabels;
    
    % find automatically labelled data
    for j=1:length(automatedLabelFolders)
       if strfind(automatedLabelFolders(j).name,currDate)
           currAutomatedFolder=automatedLabelFolders(j).name;
           break
       end
    end
    
    % load automated cluster label data
    load([automatedLabelHome '/' currAutomatedFolder '/fly' currFlyNum '/' currLobe '/' automatedLabelFilePrefix currVol '.mat'])
    
    labelComparison{i}.automated=cell(1,length(manualClusterLabels));
    labelComparison{i}.manual=cell(1,length(manualClusterLabels));
    % record automated and manual labels
    for j=1:length(clusterAssignments)
        labelComparison{i}.automated{clusterAssignments(j)}=clusterLabels{j};
    end
    for j=1:length(manualClusterLabels)
        labelComparison{i}.manual{j}=manualClusterLabels{j};
    end
    
    % create confusion matrix
    individualConfusionMatrix{i}=zeros(length(publishedOR.gh146glomerulusNames));
    % fill in confusion matrix
    for j=1:length(manualClusterLabels)
        manualLabel=[];
        automatedLabel=[];
        for k=1:length(publishedOR.gh146glomerulusNames)
            if strcmp(labelComparison{i}.manual{j},publishedOR.gh146glomerulusNames{k})
               manualLabel=k;
                break 
            end
        end
        for m=1:length(publishedOR.gh146glomerulusNames)
            if strcmp(labelComparison{i}.automated{j},publishedOR.gh146glomerulusNames{m})
               automatedLabel=m;
                break 
            end
        end
        if length(automatedLabel)>0 && length(manualLabel)>0
            individualConfusionMatrix{i}(automatedLabel,manualLabel)=1;
        end
    end
    totalConfusion=totalConfusion+individualConfusionMatrix{i};
end

%  normalize confusion matrix by number of manual labels and convert to
%  percentage
%totalConfusion=100*totalConfusion./sum(totalConfusion,1);

figure;
imagesc(totalConfusion)
xlabel('Manual Label')
ylabel('Automated Label')
set(gca,'ytick',1:length(publishedOR.gh146glomerulusNames),'yticklabel',string(publishedOR.gh146glomerulusNames),'FontSize',10)
set(gca,'xtick',1:length(publishedOR.gh146glomerulusNames),'xticklabel',string(publishedOR.gh146glomerulusNames),'FontSize',10)
xtickangle(45)
hcb=colorbar;
title(hcb,'instances')
set(gca,'FontSize',10)