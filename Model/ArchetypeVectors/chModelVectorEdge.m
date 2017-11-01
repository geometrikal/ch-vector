function [ model ] = chModelVectorEdge( ampV, angleV, N, normalise )
% CHMODELVECTOREDGE - Unweighted CH vector corresponding to multiple edges
%
% Inputs:
%
% ampV        Amplitude(s) of edge segments
% orientV     Orient(s) of edge segments
% N           Maximum RT order to use
% normalise   0* - nothing, 1 - normalise vector (optional)
%
% Multiple amplitudes and orientations can be used to create the
% vector for a multi-edge model. In this case the different model
% componentents are added together.
%
% Note that the vector is unweighted.
%
% Outputs:
%
% model       Model CH vector (unweighted)
%
%
% Written by:
%
% Ross Marchant
% James Cook University
% ross.marchant@my.jcu.edu.au
%

% Parse input arguments
if nargin < 4
    normalise = 1;
end

% Pre-allocate
model = zeros(1, 2*N+1);
edge = zeros(1, 2*N+1);

% Create single edge segment vector
for j = -N:N
    idx = j+N+1;
        edge(idx) = edge(idx) + 1 * (-1i).^(abs(j+1)) * (exp(1i*j*0));
end

% Fix 0th order
edge(N+1) = 0;

% Create vector by adding individual components
for k = 1:numel(ampV)
    model = model + ampV(k) * edge .* chRotationVector(angleV(k),N);
end

% Normalise?
if normalise == 1
    model = model ./ norm(model);
end

end


