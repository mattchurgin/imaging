clear all
fname=uigetfile;
savename=[fname(1:end-4) '_manualLabels.mat'];
load(fname)

clusterLabels=cell(1,length(clusterVolU));
for j=1:length(clusterLabels)
    clusterLabels{j}=num2str(j);
end

if exist(savename)
    disp(['cluster labels already exist for file ' fname ' . aborting'])
else
    save(savename,'clusterLabels','savename')
    disp(['initialized and saved manual cluster labels'])
end