function [ I, phase ] = chImageSinusoid( width, waveV, ampV, phaseV, orientV )
%CHIMAGESINUSOID - Image of a sinusoid

% Pre-allocate
I = zeros(width);

for j = 1:numel(waveV)    
    [x,y] = meshgrid(1:width,1:width);
    x=x-1 - floor(width/2);
    y=y-1 - floor(width/2);    
    dn = 2*pi*x./waveV(j)*cos(orientV(j)) + ...
        2*pi*y./waveV(j)*sin(orientV(j)) + phaseV(j);
    I = I + ampV(j)*cos(dn);
    phase(:,:,j) = mod(dn,2*pi);
end

