function [ f, A, theta ] = dzVector( I, N, filterType, filterParams, varargin )
%DZVECTOR - Wirtinger vector of image
%
% See chVector for inputs/outputs. The only difference is that instead of a
% Riesz transform operator, the Wirtiniger operator is used.
%
% The Wirtinger operator is essentially a complex-valued derivative of both
% axes. The Riesz transform operator is its frequency-normalised form. Thus
% instead of increasing in spatial extent with higher orders as with the
% Riesz transform, using the Wirtinger operator the spatial extent of the
% wavelets stay the same and instead the frequency band increases. Note,
% that this means tight wavelet frames are not usually possible.
%
%
% Written by:
%
% Ross Marchant
% James Cook University
% ross.marchant@my.jcu.edu.au
%

% Parse inputs
p = inputParser;
addRequired(p,'I',@ismatrix);
addRequired(p,'N',@isscalar);
addRequired(p,'filterType');
addRequired(p,'filterParams');
expPassband = {'band','high','low'};
addOptional(p,'Passband', 'band', @(x) any(validatestring(lower(x),expPassband)));
addOptional(p,'Squared', 0, @isscalar);
p.parse(I,N,filterType,filterParCHams,varargin{:})
argin = p.Results

% Ensure image is double type and single channel
I = double(I(:,:,1));

% Create FFT mesh
[ux,uy,r,th] = rtFFTMesh(size(I));

% FFT of image
Ifft = fft2(I);

% Fix to set wavelength 2 parts to zero. Not ideal as we are losing some
% information but otherwise the positive and negative orders will not be in
% conjugate.
[sr,sc] = size(I);
if mod(sr,2) == 0
    Ifft(sr/2+1,:) = 0;
end
if mod(sc,2) == 0
    Ifft(:,sc/2+1) = 0;
end

% Apply filter
filterFFT = chFilterSpectrum(size(I),filterType,filterParams,...
                             argin.Passband,argin.Squared);
Ifft = Ifft .* filterFFT;

% Wirtinger operator
RT = (ux + 1i*uy);
RT(1,1) = 0;

% 0th order of vector
f(:,:,N+1) = real(ifft2(Ifft));
A(:,:,N+1) = f(:,:,N+1);
theta(:,:,N+1) = zeros(size(f(:,:,N+1)));

% Remaining orders of vector
RTN = RT;
for j = 1:N
    % Complex form
    f(:,:,N+1+j) = ifft2(RTN.*Ifft);
    f(:,:,N+1-j) = ifft2(conj(RTN).*Ifft);
    % Polar form
    A(:,:,N+1+j) = abs(f(:,:,N+1+j));
    A(:,:,N+1-j) = abs(f(:,:,N+1-j));
    theta(:,:,N+1+j) = angle(f(:,:,N+1+j));
    theta(:,:,N+1+j) = angle(f(:,:,N+1+j));
    RTN = RTN .* RT;
end     

end

