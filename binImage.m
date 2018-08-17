function binnedImage = binImage(imageToBin,bins)
% binImage bins pixels of an image
% imageToBin is a 2-d input image
% bins is the number of pixels to bin
% output is the binned image
% by Matt Churgin

[sizex sizey]=size(imageToBin);

binnedImage=zeros(floor(sizex/bins),floor(sizey/bins));

ix=1;
for i=1:bins:(sizex-(bins-1))
    jx=1;
    for j=1:bins:(sizey-(bins-1))
        binnedImage(ix,jx)=mean(mean(imageToBin(i:(i+bins-1),j:(j+bins-1))));
        jx=jx+1;
    end
    ix=ix+1;
end