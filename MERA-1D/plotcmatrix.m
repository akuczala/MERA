function plotcmatrix(mat)
dims = size(mat);
img = zeros([dims, 3]);
%img(:,:,1) = (real(mat)+1)/2;
%img(:,:,3) = (imag(mat)+1)/2;
hue = (angle(mat)/pi+1)/2;
lum = abs(mat);

img(:,:,1) = hue.*lum;
img(:,:,3) = (1-hue).*lum;
image(img)

return