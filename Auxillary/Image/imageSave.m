function [ ] = imageSave( filename, scale, enabled )

if ~exist('enabled')
    enabled = 1;
end

if ~exist('scale')
    scale = 1;
end

if enabled ~= 0
    imageFillWindow(size(getimage), scale)
    export_fig(filename,'-native','-nocrop');
end

end

