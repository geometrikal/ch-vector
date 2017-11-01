function [  ] = imagescRange( image, crange )
%IMAGESCRANGE Summary of this function goes here
%   Detailed explanation goes here

if ~exist('crange')
    crange = [0,1];
end

imagesc(image);
caxis(crange);

end

