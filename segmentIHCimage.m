clear all
close all

filePrefix='fly15_bothALs_';
nc82name=[filePrefix 'nc82.tif'];
brpname=[filePrefix 'brpshort.tif'];
spotname=[filePrefix 'spots.tif'];
surfacename=[filePrefix 'surfaces.tif'];

% get image info
nc82info=imfinfo(nc82name);
brpinfo=imfinfo(brpname);
spotinfo=imfinfo(spotname);
surfaceinfo=imfinfo(surfacename);

nc82l=length(nc82info);
brpl=length(brpinfo);
spotl=length(spotinfo);
surfacel=length(surfaceinfo);

nc82im=zeros(nc82info(1).Width,nc82info(1).Height,nc82l);
brpim=zeros(brpinfo(1).Width,brpinfo(1).Height,brpl);
spotim=zeros(spotinfo(1).Width,spotinfo(1).Height,spotl);
surfaceim=zeros(surfaceinfo(1).Width,surfaceinfo(1).Height,surfacel);
for i=1:nc82l
    try
        nc82im(:,:,i)=imread(nc82name,'Index',i);
        
        brpim(:,:,i)=imread(brpname,'Index',i);
        
        spotim(:,:,i)=imread(spotname,'Index',i);
        
        surfaceim(:,:,i)=imread(surfacename,'Index',i);
    catch
        nc82im(:,:,i)=transpose(imread(nc82name,'Index',i));
        
        brpim(:,:,i)=transpose(imread(brpname,'Index',i));
        
        spotim(:,:,i)=transpose(imread(spotname,'Index',i));
        
        surfaceim(:,:,i)=transpose(imread(surfacename,'Index',i));
    end
    if mod(i,10)==0
        disp(['loaded ' num2str(i)])
    end
end
disp('images loaded')


gs=30;
x=[-100:1:100];
y=x;
z=-1:1; 
[xx yy zz]=meshgrid(x,y,z);
gau=1/(length(x)*gs^2)*exp(-(xx.^2+yy.^2 +zz.^2)/(2*gs^2));
disp('normalizing nc82 image')
tic
nc82norm=nc82im./(convn(nc82im,gau,'same'));

disp(['nc82 normalization done. time elapsed: ' num2str(toc/60) ' minutes'])
% 
% gs=30;
% x=[-100:1:100];
% y=x;
% z=-1:1; 
% [xx yy zz]=meshgrid(x,y,z);
% gau=1/(length(x)*gs^2)*exp(-(xx.^2+yy.^2 +zz.^2)/(2*gs^2));
% disp('normalizing brp-short image')
% tic
% brpnorm=brpim./(convn(brpim,gau,'same'));
% 
% disp(['brp normalization done. time elapsed: ' num2str(toc/60) ' minutes'])

% smooth spot image to remove outlier spots
spotimthresh=0.001;
minarea=1000000;

% disp('filling brpnorm holes')
% brpfilled=imfill(brpnorm);
% 
% disp('filling nc82norm holes')
% nc82filled=imfill(nc82norm);

% convolve binspot with gaussian
gs=30;
x=[-100:1:100];
y=x;
z=-2:2; 
[xx yy zz]=meshgrid(x,y,z);
gau=1/(length(x)*gs^2)*exp(-(xx.^2+yy.^2 +zz.^2)/(2*gs^2));
disp('smoothing spot image')
binspotsmooth=(convn(uint8(spotim),gau,'same'))>spotimthresh;
disp('pruning spot image')
binspotpruned=bwareaopen(binspotsmooth,minarea);

disp('masking nc82 image')
maskednc82=binspotpruned.*nc82norm;

bothsides=maskednc82;

nc82norm=single(nc82norm);
brpnorm=single(brpnorm);
bothsides=single(bothsides);
spotim=spotim>0;
save([filePrefix 'normed'],'nc82l','nc82norm','spotim','bothsides')

%% manually mask left and right lobes
figure
imagesc(nc82norm(:,:,30))
rmask=roipoly;
lmask=~rmask;

for i=1:nc82l
   imagesc(rmask.*nc82norm(:,:,i),[0 30])
   pause(0.05)
end
for i=1:nc82l
   imagesc(lmask.*nc82norm(:,:,i),[0 30])
   pause(0.05)
end



rightim=rmask.*maskednc82;
leftim=lmask.*maskednc82;

leftim(isnan(leftim))=0;
rightim(isnan(rightim))=0;

disp('done')

%% set threshold on slice
  
upperthresh=10; %was 10
%upperthresh=100; % was 100
lowerthresh=upperthresh; 
thresh=linspace(upperthresh,lowerthresh,size(nc82im,3));

sliceNum=10;

temp=bothsides(:,:,sliceNum);

preprocessdisksize=1; % was 3
tempn2=imopen(temp,strel('disk',preprocessdisksize));
%I = imgaussfilt(tempn2,2);
I = tempn2;
im=I<thresh(sliceNum);
im=imopen(im,strel('disk',preprocessdisksize));
im=bwareaopen(im,500); % was 100
im=~bwareaopen(~im,500); % was 100

figure
subplot(1,2,1)
imagesc(temp,[0 5*upperthresh])

% Distance transform
imb=bwdist(im);

% Blur
sigma=3;
kernel = fspecial('gaussian',4*sigma+1,sigma);
im2=imfilter(imb,kernel,'symmetric');
%im2=imb;

% Watershed transform
minSuppressionThreshold=12; 
L = watershed(imhmin(max(im2(:))-im2,minSuppressionThreshold));
lblImg = bwlabel(L&~im);

subplot(1,2,2);
%imshow(label2rgb(lblImg,'jet','k','shuffle'));
imagesc(lblImg);

%% perform watershed on each 2d slice in z-stack
  
myimtrunc=bothsides(:,:,end:-1:1);
upperthresh=9; % was 100
lowerthresh=upperthresh; 
thresh=linspace(upperthresh,lowerthresh,size(nc82im,3));

lblImg=zeros(size(rightimbrp));
disp('beginning 2d watershed')
figure
for sliceNum=1:nc82l
    
    temp=myimtrunc(:,:,sliceNum);
    
    preprocessdisksize=2; % was 3
    tempn2=imopen(temp,strel('disk',preprocessdisksize));
    %I = imgaussfilt(tempn2,2);
    I = tempn2;
    im=I<thresh(sliceNum);
    im=imopen(im,strel('disk',preprocessdisksize));
    im=bwareaopen(im,500); % was 100
    im=~bwareaopen(~im,500); % was 100
   
    % Distance transform
    imb=bwdist(im);
    
    % Blur
    sigma=3;
    kernel = fspecial('gaussian',4*sigma+1,sigma);
    im2=imfilter(imb,kernel,'symmetric');
    %im2=imb;
    
    % Watershed transform
    minSuppressionThreshold=15;
    L = watershed(imhmin(max(im2(:))-im2,minSuppressionThreshold));
    lblImg(:,:,sliceNum) = bwlabel(L&~im);
    
    subplot(1,2,1)
    imagesc(myimtrunc(:,:,sliceNum),[0 50])
    
    subplot(1,2,2)
    imagesc(lblImg(:,:,sliceNum))
    
    drawnow
    
    disp(['processed slice ' num2str(sliceNum)])
end

% perform 3d watershed transform on entire z-stack
disp('beginning 3d watershed')
upperthresh=10;
lowerthresh=upperthresh;

%myimtrunc=imgaussfilt3(rightim,[2 2 1]);

thresh=linspace(upperthresh,lowerthresh,size(myimtrunc,3)); 

disp('loading z-stack')
thresholdedim=zeros(size(myimtrunc));
for sliceNum=1:(size(myimtrunc,3))
    tempn=myimtrunc(:,:,sliceNum);
    preprocessdisksize=3; % was 3
    tempn2=imopen(tempn,strel('disk',preprocessdisksize)); 
    I = imgaussfilt3(tempn2,2);
    I = tempn2;
    im=I<thresh(sliceNum);
    im=imopen(im,strel('disk',preprocessdisksize));
    im=bwareaopen(im,500); 
    im=~bwareaopen(~im,500); 
    thresholdedim(:,:,sliceNum)=im;
end

disp('opening z-stack')
% open the image in 3d with ellipsoid
thresholdedim2=imopen(thresholdedim,strel('ball',7,3))==1;

disp('finished pre-processing z-stack')

% Distance transform
imb=bwdist(thresholdedim2);

% Blur
 %sigma=3;
 %im2=imgaussfilt3(imb,sigma);
im2=imb; % good results with no smoothing before

% Watershed transform
tic
disp('beginning watershed transform')
minSuppressionThreshold=5; % was 5.
preL=max(im2(:))-im2;
preL2=imhmin(preL,minSuppressionThreshold);
L = watershed(preL2);
lblImg3d = bwlabeln(L&~thresholdedim2);
lblImg3d(:,:,1)=0;
lblImg3d(:,:,end)=0;

figure
subplot(1,2,1)
imagesc(myimtrunc(:,:,end-5),[0 500])
subplot(1,2,2)
imagesc(lblImg3d(:,:,end-5))

disp(['finished. time elapsed: ' num2str(toc)])

save([filePrefix 'watershed.mat'],'lblImg','lblImg3d','bothsides')
%% step through labelled slices
figure
for i=1:nc82l
    subplot(2,3,1)
    imagesc(myimtrunc(:,:,i),[0 800])
    
    subplot(2,3,2)
    imagesc(thresholdedim2(:,:,i))
    
    subplot(2,3,3)
    imagesc(thresholdedim(:,:,i))
    subplot(2,3,4)
    imagesc(imb(:,:,i))
    
    subplot(2,3,6)
    imagesc(lblImg3d(:,:,i))
    
    drawnow
    %pause(0.05)
end
%% save labelled slices as gif
figure
filename='right_nc82_watershed.gif';
for i=1:size(lblImg,3)
    subplot(1,2,1)
    imagesc(myimtrunc(:,:,i),[0 100])
    subplot(1,2,2)
    imagesc(lblImg(:,:,i))
    
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')
    %text(-20,-25,['z-slice ' num2str(i)],'FontSize',20)
    
    drawnow
    
    % Capture the plot as an image
    frame = getframe;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
     % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.5);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.51);
    end
    
end

%% find voxel islands and plot

lbl2=lblImg3d;
disp('finding pixel islands')

warning('off')
areathresh=5000;

lbls=regionprops(lbl2,'all');
todelete=[];
clusterVols2=cell(1,length(lbls));
clusterInfo2=cell(1,length(lbls));
for j=1:length(lbls)
    if lbls(j).Area>areathresh 
        clusterVols2{j}=(lbl2==j);
        clusterInfo2{j}=lbls(j);
    else
        todelete=[todelete j];
    end
    if mod(j,10)==0
       disp(['processed ' num2str(j)]) 
    end
end
clusterVols2(todelete)=[];
clusterInfo2(todelete)=[];

disp(['plotting results. found ' num2str(length(clusterVols2)) ' voxel islands'])

% plot downsampled results
downsamples=4;

figure
view(3);
axis tight
camlight
lighting gouraud

hold on

rng('default')
mycmap=hsv(length(clusterVols2));
mycmap=mycmap(randperm(length(mycmap)),:);

randv=randperm(length(clusterVols2));
for i=randv
    p2=patch(isosurface(clusterVols2{i}(downsamples:downsamples:end,downsamples:downsamples:end,:)),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
    isonormals(clusterVols2{i}(downsamples:downsamples:end,downsamples:downsamples:end,:),p2)

    drawnow
end
disp('done plotting')


%% marker-controlled water shed on slices
myimtrunc=bothsides(:,:,end:-1:1);

sliceNum=50;
thresh=10;
figure
for sliceNum=70
I=myimtrunc(:,:,sliceNum);
I=imgaussfilt(I,2);

gmag = imgradient(I);
%imshow(gmag,[])

se = strel('disk',8);
Io = imopen(I,se);
% figure;
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Io)

se = strel('disk',8);
Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
% figure;
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Iobr)

Ioc = imclose(Io,se);
% figure;
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Ioc)

Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Iobrcbr)

fgm = imregionalmax(Iobrcbr>thresh);
% figure
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(fgm)

im=bwareaopen(fgm,1000);


% Distance transform
imb=imgaussfilt(bwdist(~im),2);

minSuppressionThreshold=10;

try
L = watershed(imhmin(max(imb(:))-imb,minSuppressionThreshold));
templbl = bwlabel(L&im);
% 
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm,se2);
% fgm4 = bwareaopen(fgm2,500);
% %figure;imagesc(I3)
% % 
% bw = bwareaopen(Iobrcbr>thresh,500);
% %figure;imagesc(bw)
% D = bwdist(bw);
% DL = watershed(D);
% bgm = DL == 0;
% %imagesc(bgm)
% 
% gmag2 = imimposemin(gmag, bgm | fgm4);
% 
% bw2=imgaussfilt(bwdist(gmag2>0),2);
% 
% L = watershed(imhmin(max(bw2(:))-bw2,minSuppressionThreshold));
% templbl2 = bwlabel(L&im);


subplot(1,3,1)
imagesc(I,[0 50])

subplot(1,3,2)
imagesc(templbl)

% subplot(1,3,3)
% imagesc(templbl2)
drawnow
catch
end

end

%% marker-controlled water shed on entire volume
myimtrunc=bothsides(:,:,(end-1):-1:2); % cut off first and last slices

thresh=10.5; 


tic

disp('smoothing image')
I=myimtrunc;
I=imgaussfilt3(I,[2 2 0.5]);
disp(['time elapsed: ' num2str(toc/60) ' minutes'])

% gmag = imgradient(I);
% %imshow(gmag,[])

% disp('opening image')
% %se = strel('ball',10,3);
% se = strel('ball',10,4);
% Io = imopen(I,se);
% disp(['time elapsed: ' num2str(toc/60) ' minutes'])
% figure;
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Io)

disp('opening image by reconstruction')
se = strel('ball',5,2);  % original
se = strel('ball',3,1);  % this may work better... less bleed across glomeruli
Ie = imerode(I,se);
Iobr = imreconstruct(Ie,I);
disp(['time elapsed: ' num2str(toc/60) ' minutes'])
% figure;
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Iobr)

% disp('closing image')
% se = strel('ball',3,1);
% Ioc = imclose(Io,se);
% disp(['time elapsed: ' num2str(toc/60) ' minutes'])
% figure;
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Ioc)

disp('closing by reconstruction')
se = strel('ball',10,5); % default
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
disp(['time elapsed: ' num2str(toc/60) ' minutes'])
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(Iobrcbr)

disp('finding regional max')
fgm = imregionalmax(Iobrcbr>thresh); % default
disp(['time elapsed: ' num2str(toc/60) ' minutes'])
% figure
% subplot(1,2,1);imagesc(I)
% subplot(1,2,2);imagesc(fgm)

im=fgm;
%im=bwareaopen(im,1000);
%im=~bwareaopen(~im,5000);


% Distance transform
disp('calculating distance transform')

imb=imgaussfilt3(bwdist(~im),[2 2 0.5]);

disp(['time elapsed: ' num2str(toc/60) ' minutes'])

minSuppressionThreshold=4; % was 5

disp('calculating watershed transform')
L = watershed(imhmin(max(imb(:))-imb,minSuppressionThreshold));
templbl = bwlabeln(L&im);
disp(['time elapsed: ' num2str(toc/60) ' minutes'])

disp('finished watershed')

% find pixel islands and plot
lbl2=templbl;
disp('finding pixel islands')

warning('off')
areathresh=5000;

lbls=regionprops(lbl2,'all');
todelete=[];
clusterVols2=cell(1,length(lbls));
clusterInfo2=cell(1,length(lbls));
for j=1:length(lbls)
    if lbls(j).Area>areathresh 
        clusterVols2{j}=(lbl2==j);
        clusterInfo2{j}=lbls(j);
    else
        todelete=[todelete j];
    end
    if mod(j,10)==0
       disp(['processed ' num2str(j)]) 
    end
end
clusterVols2(todelete)=[];
clusterInfo2(todelete)=[];

disp(['plotting results. found ' num2str(length(clusterVols2)) ' voxel islands'])

% plot downsampled results
downsamples=4;

figure
view(3);
axis tight
camlight
lighting gouraud

hold on

rng('default')
mycmap=hsv(length(clusterVols2));
mycmap=mycmap(randperm(length(mycmap)),:);

randv=randperm(length(clusterVols2));
j=1;
for i=randv
    p2=patch(isosurface(clusterVols2{i}(downsamples:downsamples:end,downsamples:downsamples:end,:)),'FaceColor',mycmap(i,:),'EdgeColor','none','FaceAlpha',0.8);
    isonormals(clusterVols2{i}(downsamples:downsamples:end,downsamples:downsamples:end,:),p2)

    drawnow
    disp(['plotted island ' num2str(j) ' of ' num2str(length(clusterVols2))])
    j=j+1;
end
disp('done plotting')

save([filePrefix 'watershednew2.mat'],'templbl','bothsides','thresh','minSuppressionThreshold')
disp(['done. total time elapsed: ' num2str(toc/60) ' minutes.'])



%% step through labelled slices
% 
% figure
% for i=1:size(templbl,3)
% 
%     imagesc(nc82norm(:,:,i),[0 30])
%     pause
% end
figure
for i=1:size(templbl,3)
    %subplot(1,2,1)
    %imagesc(myimtrunc(:,:,i),[0 30])
    
%     subplot(1,3,2)
%     imagesc(lblImg3d(:,:,i))
%     
    %subplot(1,2,2)
    imagesc(templbl(:,:,i))

    

    drawnow
   % pause(0.05)
   pause
end