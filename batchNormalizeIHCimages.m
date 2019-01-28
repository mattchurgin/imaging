% batch normalize ihc images
clear all
close all

folders{1}='190117_orcoOXZbrpshort';
folders{2}='190122_orcoOXZ_brpshort';
homeDir='D:\HCBI';

for i=1:length(folders)
    cd(folders{i})
    currfolders=dir(pwd);
    currfolders=currfolders(3:end);
    for j=1:length(currfolders)
        cd(currfolders(j).name)
        
        filePrefix=[currfolders(j).name '_bothALs_'];
        nc82name=[filePrefix 'nc82.tif'];
        spotname=[filePrefix 'spots.tif'];
        surfacename=[filePrefix 'surfaces.tif'];
        
        % get image info
        nc82info=imfinfo(nc82name);
        spotinfo=imfinfo(spotname);
        surfaceinfo=imfinfo(surfacename);
        
        nc82l=length(nc82info);
        spotl=length(spotinfo);
        surfacel=length(surfaceinfo);
        
        nc82im=zeros(nc82info(1).Width,nc82info(1).Height,nc82l);
        spotim=zeros(spotinfo(1).Width,spotinfo(1).Height,spotl);
        surfaceim=zeros(surfaceinfo(1).Width,surfaceinfo(1).Height,surfacel);
        for i=1:nc82l
            try
                nc82im(:,:,i)=imread(nc82name,'Index',i);
                
                spotim(:,:,i)=imread(spotname,'Index',i);
                
                surfaceim(:,:,i)=imread(surfacename,'Index',i);
            catch
                nc82im(:,:,i)=transpose(imread(nc82name,'Index',i));
                
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
        disp(['normalization convolution done. time elapsed: ' num2str(toc/60) ' minutes'])

        save([currfolders(j).name '_normed'],'nc82norm','spotim','surfaceim','nc82info')
        cd ..
    end
    cd(homeDir)
end