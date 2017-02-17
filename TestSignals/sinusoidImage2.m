function [ I, p ] = sinusoidImage2( height, width, waveV, ampV, phaseV, orientV )
%SINUSOIDIMAGE Summary of this function goes here
%   Detailed explanation goes here

I = zeros([height,width]);


for j = 1:numel(waveV)
    
    [y,x] = meshgrid(1:width,1:height);
    x=x-1 - floor(height/2);
    y=y-1 - floor(width/2);
    
    dn = 2*pi*x./waveV(j)*cos(orientV(j)) + ...
        2*pi*y./waveV(j)*sin(orientV(j)) + phaseV(j);
    I = I + ampV(j)*cos(dn);
    p(:,:,j) = mod(dn,2*pi);

end

