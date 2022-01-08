function nImg = normalizeImageToDAPI(img,dapi,rSum)

newimg = im2double(img);
dapiimg = im2double(dapi);

% flt = fspecial('disk',rSum);
flt = fspecial('gaussian',rSum,5) ;
flt(flt > 0)= 1.0;

newimg = imfilter(newimg,flt,'replicate');
dapiimg = imfilter(dapiimg,flt,'replicate');

newimg = newimg./dapiimg;

newimg(dapiimg == 0) = 0;

nImg = newimg;