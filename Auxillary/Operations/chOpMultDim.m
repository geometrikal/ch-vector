function [ result ] = chOpMultDim( matrix, vector, dimension )
%CHOPMULTDIM - Multiply a matrix by a vector along a certain dimension
%
% Inputs:
%
% matrix        Matrix, or cell array of matrices, to multiply
% vector        Vector to multiply
% dimension     Dimension for multiplication (optional)
%
% If no dimension is specified, the multiplication will occur along the
% highest dimension whose size matches the size of the vector
%
% Outputs:
% 
% result        Result of multiplication
%

% Parse input arguments
if nargin < 3
    if ~iscell(matrix)
        s = size(matrix);
    else
        s = size(matrix{1});
    end
    n = numel(vector);
    idx = find(s == n);
    if numel(idx) > 0
        dimension = idx(end);
    else
        error('None of the dimensions of the matrix match that of the vector');
    end
end

vector = vector(:);
shape = ones(1,numel(size(matrix)));
shape(dimension) = numel(vector);

vector = reshape(vector,shape);

if ~iscell(matrix)
    result = bsxfun(@times, matrix, vector);
else
    for k = 1:numel(matrix)
        result{k} = bsxfun(@times, matrix{k}, vector);
    end
end

end

