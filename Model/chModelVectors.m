function [ U ] = chModelVectors( modelType, N, weights, params )
%CHMODELVECTORS - Creates set of CH vectors in matrix representing a model
%
% Inputs:
%
% modelType         The type of model
% - sinusoidal      Sinusoidal model (odd and even)
% - line            Line segment
% - edge            Edge segment
% - wedge           Wedge segment (not implemented yet)
%                   params is angle in radians
% - halfSinusoidal  Half-sinusoidal model (line and edge)
%
% Outputs:
%
% U                 Matrix of model CH vectors. The columns correspond to 
%                   vectors, and the rows to CH orders.
%
% Written by:
%
% Ross Marchant
% James Cook University
% ross.marchant@my.jcu.edu.au
%

switch lower(modelType)
    
    case {'sinu','sinusoidal','sinusoid'}
        U(:,1) = chModelVectorSinusoid(1,0,0,N,1);
        U(:,2) = chModelVectorSinusoid(1,pi/2,0,N,1);
        
    case {'line'}
        U(:,1) = chModelVectorLine(1,0,N,1);
        
    case {'edge'}
        U(:,1) = chModelVectorEdge(1,0,N,1);
        
    case {'half','halfsinusoidal','halfsinusoid'}
        U(:,1) = chModelVectorLine(1,0,N,1);
        U(:,2) = chModelVectorEdge(1,0,N,1);
        
end

% Check weights
if size(weights,1) > size(weights,2)
    weights = permute(weights,[2,1]);
end

% Apply weighting
U = bsxfun(@times, U, weights.');
U = bsxfun(@rdivide,U,chOpNormDim(U,1));

end

