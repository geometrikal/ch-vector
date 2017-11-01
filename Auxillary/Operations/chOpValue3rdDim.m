function [ m ] = chOpValue3rdDim( matrix, idx )
%CHOPVALUE3RDDIM 3D to 2D matrix, selecting one value on 3rd dim.

[ri,ci] = ndgrid(1:size(matrix,1),1:size(matrix,2));
ri = repmat(ri,1,1,size(idx,3));
ci = repmat(ci,1,1,size(idx,3));
m = matrix(sub2ind(size(matrix),ri,ci,idx));

end

