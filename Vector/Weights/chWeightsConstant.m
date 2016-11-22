function [ c ] = chWeightsConstant( N )
%CHWEIGHTSCONSTANT - Constant for normalisation of width parameter with N
% when calculating weights

c = 5.64/N - 6.57/N^2;

end

