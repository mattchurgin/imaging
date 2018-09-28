function [] = plotALResponse(grnResponse,t,odors)
% plots calcium responses contained in grnResponse
% t specifies the time vector
% odors is a vector specifying which odors to present the response for
% if odors is not specified, all odor responses will be shown

if nargin<3
    odors=1:size(grnResponse,1);
end

figure
for i=1:length(odors)
    subplot(1,length(odors),i)
    
    imagesc(t,1:size(grnResponse,2),squeeze(grnResponse(odors(i),:,:)),[0 1])
    title(['odor ' num2str(i)])
    if i==1
        xlabel('Time (s)')
        ylabel('Cluster #')
    end
    
    title(['Odor #' num2str(i-1)])
    set(gca,'YTick','')
    
    set(gca,'FontSize',15)
end
