function [ux,uy,r,th] = chFFTMesh(sizeV)
%CHFFTMESH - Create DFT cartesian and polar coordinates
%
% Inputs:
% 
% sizeV     [width,height] of spectrum (same as image size)
%
% Outputs:
%
% ux        x in cartesian coords
% uy        y in cartesian coords
% r         radius in polar coords  
% th        angle (rad CCW) in polar coords

rows = sizeV(1);
cols = sizeV(2);


[ux, uy] = meshgrid(...
    ([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
    ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2))...
    );
ux = ifftshift(ux);   % Quadrant shift to put 0 frequency at the corners
uy = ifftshift(uy);

% Convert to polar coordinates
th = atan2(uy,ux);
r = sqrt(ux.^2 + uy.^2);

end
