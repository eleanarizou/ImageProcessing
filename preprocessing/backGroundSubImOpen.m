function img_bgsub = backGroundSubImOpen(img,BGradius)

if ~exist('BGradius','var')
    BGradius = 20;
end

img_bgsub = zeros(size(img));

for ii = 1:size(img_bgsub,3)
    bg = imopen(img(:,:,ii),strel('disk',BGradius));
    img_bgsub(:,:,ii) = imsubtract(img(:,:,ii),bg);
end