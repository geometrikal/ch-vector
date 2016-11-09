function [ f, A, theta ] = chVectorMultiScale( I, N, numScales, filterType, filterParams, varargin )
%CHVECTORMULTISCALE - CH vector of image at multiple scales
%
% Inputs:
%
% I             Image
% N             Maximum circular harmonic order
% numScales     Number of scales 
% filterType    filter type       - see chFilterSpectrum
% filterParams  filter parameters - see chFilterSpectrum
%
% Optional inputs:
%
% Factor        2* - the factor to multiply the wavelength by at each scale
% Squared       0* - nothing, 1 - square the filter freq. response
% Passband      band*, low, high
%
%               The first scale will be at the wavelength specified in the
%               filter parameters, for subsequent scales the wavelength
%               will be multiplier by the factor each time. The default
%               factor is 2 which gives filters spaced an octave apart.
%
% Outputs:
% 
% f             Cells of CH vector matrix where k is the scale order:
%               f{k}(:,:,1)     : -Nth order
%               f{k}(:,:,N+1)   :  0th order
%               f{k}(:,:,2*N+1) :  Nth order
%
% A             Amplitude of polar form of each scale and order - abs(f{k}(:,:,n))
% theta         Orientation of polar form of each scale and order - angle(f{k}(:,:,n))
%
% Example use:  [f, A, theta] = chVectorMultiScale(I,9,4,'lognormal',[4,0.65])
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
addOptional(p,'Factor', 2, @isscalar);
p.parse(I,N,filterType,filterParams,varargin{:})
argin = p.Results

% Ensure image is double type and single channel
I = double(I(:,:,1));

% Create FFT mesh
[ux,uy,r,th] = chFFTMesh(size(I));

% FFT of image
IfftOrig = fft2(I);

% Fix to set wavelength 2 parts to zero. Not ideal as we are losing some
% information but otherwise the positive and negative orders will not be in
% conjugate.
[sr,sc] = size(I);
if mod(sr,2) == 0
    IfftOrig(sr/2+1,:) = 0;
end
if mod(sc,2) == 0
    IfftOrig(:,sc/2+1) = 0;
end

% CH operator
RT = (ux + 1i*uy) ./ abs(ux + 1i*uy);
RT(1,1) = 0;

for k = 1:numScales
    % Apply filter
    factor = [argin.Factor, k-1];
    filterFFT = chFilterSpectrum(size(I),filterType,filterParams,...
                                 argin.Passband,argin.Squared,factor);
    Ifft = IfftOrig .* filterFFT;

    % 0th order of CH vector
    f{k}(:,:,N+1) = real(ifft2(Ifft));
    A{k}(:,:,N+1) = f{k}(:,:,N+1);
    theta{k}(:,:,N+1) = zeros(size(f{k}(:,:,N+1)));

    % Remaining orders of CH vector
    RTN = RT;
    for j = 1:N
        % Complex form
        f{k}(:,:,N+1+j) = ifft2(RTN.*Ifft);
        f{k}(:,:,N+1-j) = ifft2(conj(RTN).*Ifft);
        % Polar form
        A{k}(:,:,N+1+j) = abs(f{k}(:,:,N+1+j));
        A{k}(:,:,N+1-j) = abs(f{k}(:,:,N+1-j));
        theta{k}(:,:,N+1+j) = angle(f{k}(:,:,N+1+j));
        theta{k}(:,:,N+1+j) = angle(f{k}(:,:,N+1+j));
        RTN = RTN .* RT;
    end    
end

end

