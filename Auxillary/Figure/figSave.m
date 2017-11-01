function [ ] = figSave( filename, dpi, enabled )
%FIGSAVE Summary of this function goes here
%   Detailed explanation goes here

if ~exist('dpi')
    dpi = 300;
end

if ~exist('enabled')
    enabled = 1;
end

if enabled ~= 0
    set(gcf,'color','w');
    export_fig(filename,'-transparent',['-r' num2str(dpi)],'-nocrop');
end

end

