function out_image = GE_image_rotate(in_image)

%PJN - takes an image created from GE raw data and rotates it so that axial
%is the first dimension, is correctly oriented (left on left), and goes
%from top to bottom

RotImage = zeros(size(in_image));

for i = 1:size(in_image,1)
    RotImage(i,:,:) = fliplr(rot90(squeeze(in_image(i,:,:))));
end

RotImage = flip(RotImage,1);

out_image = RotImage; 