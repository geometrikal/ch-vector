function [ diff ] = chOpAngleDifference( x, y, use_polarity )
%RTANGLEDIFFERENCE Summary of this function goes here
%   Detailed explanation goes here

diff = abs(atan2(sin(x-y), cos(x-y)));

if ~exist('use_polarity')
    use_polarity = 0;
end

if use_polarity
    diff = atan2(sin(x-y), cos(x-y));
%     polarity = ones(size(x))
%     polarity((y < x) | ( y >= x + pi)) = -1;
%     diff = diff .* polarity;
end

end

