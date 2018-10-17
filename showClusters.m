function [] = showClusters(clusterVolU,clusterInfoU,labels)
% showClusters plots clusters contained in the cell array clusterVolU in 3d
% labels is a cell array of cluster labels.  if labels is not passed in,
% clusters will be labelled with numbers

if nargin<3
   labels=[]; 
end
figure
view(3);
axis tight
camlight
lighting gouraud

hold on

rng('default')
mycmap=hsv(length(clusterVolU));
mycmap=mycmap(randperm(length(clusterVolU)),:);

for i=1:length(clusterVolU)
    p2=patch(isosurface(clusterVolU{i}),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.4);
    isonormals(clusterVolU{i},p2)
    if length(labels)<1
        text(clusterInfoU{i}.Centroid(1),clusterInfoU{i}.Centroid(2),clusterInfoU{i}.Centroid(3),num2str(i),'FontSize',20,'FontWeight','Bold')
    else
        text(clusterInfoU{i}.Centroid(1),clusterInfoU{i}.Centroid(2),clusterInfoU{i}.Centroid(3),labels{i},'FontSize',20,'FontWeight','Bold')
    end
    %drawnow
    %pause
end