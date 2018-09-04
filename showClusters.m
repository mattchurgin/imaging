function [] = showClusters(clusterVolU,clusterInfoU)
% showClusters plots clusters contained in the cell array clusterVolU in 3d

figure
view(3);
axis tight
camlight
lighting gouraud

hold on

for i=1:length(clusterVolU)
    p2=patch(isosurface(clusterVolU{i}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(clusterVolU{i},p2)
    text(clusterInfoU{i}.Centroid(1),clusterInfoU{i}.Centroid(2),clusterInfoU{i}.Centroid(3),num2str(i),'FontSize',15)
end