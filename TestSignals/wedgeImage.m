function [ I, rc, cc ] = wedgeImage( r, c, ampV, angleV, occluded, smooth, dbg )
%lineImage Creates image with radiating lines
%
%   Inputs:
%
%   r           rows
%   c           columns
%   ampV        vector of line amplitudes
%   angleV      vector of line angles in degrees
%   occluded    1 = image is max value of component lines
%               0 = image is addition of component lines
%   smooth      gaussian filter smooth amount
%   dbg         optional - enable debugging output
%
%   Outputs:
%
%   I           image
%   rc          row of line centre
%   cc          column of line centre

if nargin < 6
    dbg = 0;
end

% r = r * 4;
% c = c * 4;
% 
I = zeros(r,c);
rc = floor(r/2);
cc = floor(c/2);
rfix = mod(r,2);
cfix = mod(c,2);

angleV = mod(angleV / 180 *pi, 2*pi);

[rg,cg] = ndgrid(-rc:rc-1+rfix,-cc:cc-1+cfix);
th = atan2(rg,cg);
th = mod(th,2*pi);
%imagesc(th); pause

for k = 1:numel(angleV)/2
   
    if angleV(2*k) > angleV(2*k-1)
        temp = (th > angleV(2*k-1)) .* (th < angleV(2*k));% .*ampV(k);
    else
        temp = (th > angleV(2*k-1)) + (th < angleV(2*k));% .*ampV(k);
    end
    
%     temp = temp ./ sum(temp(:));
%     
%     temp2= zeros(r,c);
%     temp2 = temp;
%     temp2(1:end-1,:) = temp2(1:end-1,:) | temp(2:end,:);
%     temp2(:,1:end-1) = temp2(:,1:end-1) | temp(:,2:end);
%     temp2(2:end,:) = temp(1:end-1,:) | temp2(2:end,:);
%     temp2(:,2:end) = temp(:,1:end-1) | temp2(:,2:end);
%     temp =temp2;
    
    
    temp = temp .* ampV(k);
    
    if occluded == 1
        I(temp > I) = temp(temp > I);
    else
        I = I + temp;
    end
    
    if dbg == 1
        imagesc(I); colormap gray; pause;
    end

end

if smooth > 0
    %h = fspecial('gaussian',smooth*4,smooth*4/8);
    h = fspecial('gaussian',smooth*3,smooth/2);
    I = imfilter(I,h);
    %temp = temp ./ max(temp(:));
end

% I = imresize(I,0.25);
% rc =rc /4;
% cc = cc/4;


