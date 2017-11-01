function [ ] = imageFillWindow(imageSize, scale)
% IMAGEFILLWINDOW Fills the current figure window with the image (or
% otherwise) in the axes.

if nargin < 1
    imageSize = size(getimage);
end

if nargin < 2
    scale = 1;
end

daspect([1 1 1])
set(gcf,'Units','pixels','Position',[0 0 imageSize(2)*scale imageSize(1)*scale]);
set(gca,'Units','normalized','Position',[0 0 1 1]);
axis(gca,'off');

end

