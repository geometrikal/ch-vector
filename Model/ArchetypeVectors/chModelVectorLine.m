function [ model ] = chModelVectorLine( ampV, angleV, N, normalise )
% CHMODELVECTORLINE - Unweighted CH vector corresponding to multiple lines
%
% Inputs:
%
% ampV        Amplitude(s) of line segments
% orientV     Orient(s) of line segments
% N           Maximum RT order to use
% normalise   0* - nothing, 1 - normalise vector (optional)
%
% Multiple amplitudes and orientations can be used to create the
% vector for a multi-line model. In this case the different model
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

% Create vector
for j = -N:N
    idx = j+N+1;
    for k = 1:numel(ampV)
        model(idx) = model(idx) + ampV(k) * ...
            (-1*1i)^abs(j) * exp(1i*j*angleV(k));
    end
end

% Normalise?
if normalise == 1
    model = model ./ norm(model);
end

end


