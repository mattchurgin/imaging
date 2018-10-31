function [clusterVolNew clusterInfoNew] = mergeLikeClusters(clusterVolU,clusterInfoU,grnResponse)

physDistThreshold=30;
corrThreshold=0.8;

clusterCorrelation=NaN*zeros(size(grnResponse,2),size(grnResponse,2));
physDist=NaN*zeros(size(grnResponse,2),size(grnResponse,2));
highCorr=zeros(size(grnResponse,2),size(grnResponse,2));
lowDist=zeros(size(grnResponse,2),size(grnResponse,2));
mergeMatrix=zeros(size(grnResponse,2),size(grnResponse,2));
for i=1:size(grnResponse,2)
    temp1=squeeze(grnResponse(:,i,:));
    for j=1:size(grnResponse,2)
        if i~=j
            temp2=squeeze(grnResponse(:,j,:));
            tempcorr=corr(temp1(:),temp2(:),'Type','Spearman');
            clusterCorrelation(i,j)=tempcorr;
            physDist(i,j)=sqrt(sum((clusterInfoU{i}.Centroid-clusterInfoU{j}.Centroid).^2));
            
            highCorr(i,j)=clusterCorrelation(i,j)>corrThreshold;
            lowDist(i,j)=physDist(i,j)<physDistThreshold;
        end
    end
end
%
% figure;
% imagesc(clusterCorrelation)
% figure;
% imagesc(clusterCorrelation./(1+physDist),[0 0.01])


mergeMatrix=highCorr.*lowDist;
%
% figure;imagesc(highCorr)
% figure;imagesc(lowDist)
% figure;imagesc(highCorr.*lowDist)

for i=1:size(grnResponse,2)
    temp=find(mergeMatrix(i,:));
    
    for j=1:length(temp)
        for k=1:length(temp)
            if j~=k
                mergeMatrix(temp(j),temp(k))=1;
                mergeMatrix(temp(k),temp(j))=1;
            end
        end
    end
    temp=find(mergeMatrix(i,:));
end
for i=1:size(grnResponse,2)
    ms{i}=find(mergeMatrix(i,:));
end
%figure;imagesc(mergeMatrix)

toNotMerge=find(cellfun(@isempty,ms));

mergebins=cell(1,1);
nmerges=1;
for i=1:length(ms)
    if ms{i}
        mergebins{nmerges}=[i ms{i}];
        
        % don't repeat merge
        for j=1:length(ms{i})
            ms{ms{i}(j)}=[];
        end
        nmerges=nmerges+1;
    end
end

if nmerges>1
    clusterVolNew=cell(1,length(toNotMerge)+length(mergebins));
    clusterInfoNew=cell(1,length(toNotMerge)+length(mergebins));
    for i=1:length(toNotMerge)
        clusterVolNew{i}=clusterVolU{toNotMerge(i)};
        clusterInfoNew{i}=clusterInfoU{toNotMerge(i)};
    end
    for i=1:length(mergebins)
        tempvol=clusterVolU{mergebins{i}(1)};
        for j=2:length(mergebins{i})
            tempvol=tempvol+clusterVolU{mergebins{i}(j)};
        end
        tempvol=tempvol>0;
        
        CC=bwconncomp(tempvol);
        if CC.NumObjects>1
            for toremove=2:CC.NumObjects
                tempvol(CC.PixelIdxList{toremove})=0;
            end
        end
        
        clusterVolNew{length(toNotMerge)+i}=tempvol;
        clusterInfoNew{length(toNotMerge)+i}=regionprops(tempvol,'all');
    end
else
    clusterVolNew=clusterVolU;
    clusterInfoNew=clusterInfoU;
end