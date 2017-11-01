function [ f, A, theta ] = chVector( I, N, filterType, filterParams, varargin )
%CHVECTOR - CH vector of image
%
% Inputs:
%
% I             Image
% N             Maximum circular harmonic order
% filterType    filter type       - see chFilterSpectrum
% filterParams  filter parameters - see chFilterSpectrum
%
% Optional inputs:
%
% Squared       0* - nothing, 1 - square the filter freq. response
% Passband      band*, low, high
%
% Outputs:
% 
% f             CH vector matrix
%               f(:,:,1)     : -Nth order
%               f(:,:,N+1)   :  0th order
%               f(:,:,2*N+1) :  Nth order
%
% A             Amplitude of polar form of each order - abs(f(:,:,n))
% theta         Orientation of polar form of each order - angle(f(:,:,n))
%
% Example use:  [f, A, theta] = chVector(I,9,'lognormal',[32,0.65])
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
p.parse(I,N,filterType,filterParams,varargin{:})
argin = p.Results;

% Ensure image is double type and single channel
I = double(I(:,:,1));

% Create FFT mesh
[ux,uy,r,th] = chFFTMesh(size(I));

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

% CH operator
RT = (ux + 1i*uy) ./ abs(ux + 1i*uy);
RT(1,1) = 0;

% 0th order of CH vector
f(:,:,N+1) = real(ifft2(Ifft));
A(:,:,N+1) = f(:,:,N+1);
theta(:,:,N+1) = zeros(size(f(:,:,N+1)));

% Remaining orders of CH vector
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

