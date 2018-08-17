function [kmeansOut, c, sumd, d] = neuronPCA_3d(inputUnwrappedImage,imageSize,k)

kmeansOut= kmeans(inputUnwrappedImage,k,'MaxIter',500,'Display','iter','Distance','sqeuclidean');

% restore original image shape
kmeansOut = reshape(kmeansOut,[imageSize(1),imageSize(2),imageSize(3)]);

