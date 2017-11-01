function [ clr ] = figMatlabColours( idx )
%FIGMATLABCOLOURS Summary of this function goes here
%   Detailed explanation goes here

clr =  get(groot,'DefaultAxesColorOrder');
clr = clr(idx,:);

end

