function [ weights ] = chWeightsSinusoid( N, width, mode )
%CHWEIGHTSSINUSOID - Weighting matrix for sinusoidal model that ensures odd
%                    and even components are weighted equally.
%
% Inputs:
%
% N         Maximum CH order to use
% width     Width in either radians (mode = 0) or normalised to N (mode = 1*)
%           For normalised (mode = 1) the following widths result in:
%           0 - each odd order (and each order) are weighted the same. the
%           angular respons has larger lobes
%           1 - the higher orders are weighted less, resulting in smaller
%           angular lobes
%           2 - the angular response is very smooth
%           (Note: the width is a continuous parameter, e.g. width = 1.3 is
%           perfectly fine)
% mode      0 - absolute width, 1* - normalised
%
% Start with width = 0 and adjust when necessary. A width of 0 gives the
% best noise response for a single sinusoidal model, but in multiple
% component models a higher width is sometimes better so that off-axis
% components have less influence.
%
% Outputs:
%weights
% weights   The weighting vector (1 x 2N+1)
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
    mode = 1;
end

% Use normalisation?
if mode == 1
    width = width * chWeightsConstant(N);
end

% Width cannot be 0
if width == 0
    width = 0.00001;
end

N2 = 2*N+1;

% Matrix indices
n1 = repmat(-N:N,[2*N+1,1]);
n2 = repmat([-N:N]',[1,2*N+1]);
dn = n1 - n2;

% Solve for window
V = 2*sin(width*dn/2) ./ dn;
V(isnan(V)) = 0;
idx = eye(N2);
V = V + 2*width * eye(N2);
V(mod(dn,2) == 1) = 0;

% Calculate eigenvectors
[eVectors,~] = eig(V);

% Last two are the odd and even weighted sinusoidal CH vectors. Use their
% absolute values to get the weights.
weights = sum(abs(eVectors(:,end-1:end)) / sqrt(2),2);
weights = weights.';

end

