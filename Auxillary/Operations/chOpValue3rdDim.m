function [ m ] = rtValue3rdDim( matrix, idx )
%RTVALUE3RDDIM Summary of this function goes here
%   Detailed explanation goes here

[ri,ci] = ndgrid(1:size(matrix,1),1:size(matrix,2));
ri = repmat(ri,1,1,size(idx,3));
ci = repmat(ci,1,1,size(idx,3));
m = matrix(sub2ind(size(matrix),ri,ci,idx));

end

