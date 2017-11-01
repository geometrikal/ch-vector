function [ R ] = chRotationVector( theta, N, dimension)
%CHROTATATIONVECTOR - Vector for rotating CH vectors
%
% Inputs:
%
% theta       Rotation amount in radians
% N           Maximum RT order to use
% dimension   Dimension along which to rotate (default = 2)
%
% Outputs:
%
% R           Vector
%
% Usage:
% 
% Single vector, f (1 x 2N+1) 
% fRotated = f .*  chRotationVector(theta, N);
%
% Matrix of image CH vectors, f (rows x cols x 2N+1)
% fRotated = bsxfun(@times, f, chRotationVector(theta, N, 3));
%
%
% Written by:
%
% Ross Marchant
% James Cook University
% ross.marchant@my.jcu.edu.au
%

% Parse input arguments
if nargin < 3
    dimension = 2;
end

% Shape to correct dimension
n = (-N:N);
shape = ones(1,4);
shape(dimension) = numel(n);
n = reshape(n,shape);

% Create vector
R = exp(bsxfun(@times,n,1i*theta));

end

