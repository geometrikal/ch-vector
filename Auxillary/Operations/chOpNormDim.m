function [ result ] = chOpNormDim( matrix, dimension )
%CHOPNORMDIM - Calculate norm along dimension
%
% Inputs:
% 
% matrix        matrix to calculate norm
% dimension     dimension along which to calculate
%
% Outputs:
%
% result        Norm along dimension
%

result = real(sqrt(sum(matrix.*conj(matrix),dimension)));

end

