function [ model ] = chModelVectorSinusoid( ampV, phaseV, orientV, N, normalise )
% CHMODELVECTORSINUSOID - Unweighted CH vector corresponding to sinusoidal model
%
% Inputs:
%
% ampV        Amplitude(s)
% phaseV      Phase(s)
% orientV     Orient(s)
% N           Maximum RT order to use
% normalise   0* - nothing, 1 - normalise vector (optional)
%
% Multiple amplitudes, phases, and orientations can be used to create the
% vector for a multi-sinusoidal model. In this case the different model
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
        % Even order
        if mod(j,2) == 0
            model(idx) = model(idx) + ...
                ampV(k) * exp(1i*j*orientV(k)) * cos(phaseV(k));
        % Odd order
        else
            model(idx) = model(idx) + ...
                ampV(k) * exp(1i*j*orientV(k)) * 1i*sin(phaseV(k));
        end
    end
end

% Normalise?
if normalise == 1
    model = model ./ norm(model);
end

end


