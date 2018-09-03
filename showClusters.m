function [] = showClusters(clusterVolU)
% showClusters plots clusters contained in the cell array clusterVolU in 3d
view(3);
axis tight
camlight
lighting gouraud
%lighting phong
hold on

for i=1:length(clusterVolU)
    p2=patch(isosurface(clusterVolU{i}),'FaceColor',rand(1,3),'EdgeColor','none','FaceAlpha',0.3);
    isonormals(clusterVolU{i},p2)
end