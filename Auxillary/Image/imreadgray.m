function [ I ] = imreadgray( filename )
%CHLOADIMAGEGRAYSCALE - Load image, convert to grayscale and double format

I = imread(filename);

if ndims(I) > 2
    I = rgb2gray(I);
end

I = double(I);

end

