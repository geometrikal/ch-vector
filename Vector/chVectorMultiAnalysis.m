function [ f ] = chVectorMultiAnalysis( I, N, numScales, filterType, filterParams, varargin )
%CHVECTORMULTIANALYSIS - Multi-scale CH vector wavelet decompostion of image
%
% Inputs:
%
% I             Image
% N             Maximum circular harmonic order
% filterType    Filter type       - see chFilterSpectrum
%               Filter should be a wavelet for tight frame property
% filterParams  Filter parameters - see chFilterSpectrum
%
% Optional inputs:
%
% Weights       CH vector weights [w_-N, ..., w_N] - [1,...,1] default
% Squared       0* - nothing, 1 - square the filter freq. response
% Factor        See chFilterSpectrum
% Subsample     0* - none, 1 - subsample by 2 after each scale. 
%               Requires image dimensions divisible by 2^numScales
% Debug         0* - nothing, 1 - show debugging information
%
% Outputs:
% 
% f             CH vector matrix for each scale
%               f{k}(:,:,1)     : -Nth order
%               f{k}(:,:,N+1)   :  0th order
%               f{k}(:,:,2*N+1) :  Nth order
%               where k is the scale
%               f{N+1} is NOT a matrix however, it is the low pass response
%
%
% Example use:  f = chVectorMultiAnalysis(I,7,4,'vow',4)
%               will split the image into 5 scales: 1 highpass with centre
%               wavelength of 4, 3 bandpass, and 1 lowpass. The CH vector
%               will be calculated on all except the lowpass. 
%
%               Note that the wavelength (first value of filterParams)
%               should generally be 4 so that the first scale includes a
%               decent portion of the spectrum.
%
%               Weights should be set so that the sum of squared values add
%               to 1 for a Parseval tight frame.
%
%
% Written by:
%
% Ross Marchant
% James Cook University
% ross.marchant@my.jcu.edu.au
%
% TO DO:
%
% Change f to f.ch, f.lp, f.hp

% Parse inputs
p = inputParser;
addRequired(p,'I',@ismatrix);
addRequired(p,'N',@isscalar);
addRequired(p,'numScales');
addRequired(p,'filterType');
addRequired(p,'filterParams');
addOptional(p,'Squared', 0, @isscalar);
addOptional(p,'Weights', ones(1,2*N+1), @(x) (numel(x) == 2*N+1));
addOptional(p,'Factor', 2, @isscalar);
addOptional(p,'Subsample', 0, @isscalar);
addOptional(p,'Debug', 0, @isscalar);
p.parse(I,N,numScales,filterType,filterParams,varargin{:})
argin = p.Results
db = argin.Debug;

% Ensure image is double type and single channel
I = double(I(:,:,1));

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

% Process each scale
for ns = 1:numScales
    
    % *********************************
    % Make the Riesz tranform filter
    % *********************************
    if ns == 1 || argin.Subsample > 0
        [ux,uy] = chFFTMesh(size(I));
        RT = (ux + 1i*uy) ./ abs(ux + 1i*uy);
        RT(1,1) = 0;
    end
    
    % *********************************
    % Make the detail filter
    % *********************************
    % High pass if initial OR subsampling
    % Band pass if not subsampling or just one scale
    if ns == 1 || argin.Subsample > 0
        passband = 'high'
    else
        passband = 'band';
    end
    if numScales == 1
        passband = 'band'
    end
    % Choose factor - see chFilterSpectrum
    % factor(1) = division of spectrum
    % factor(2) = which bandpass filter in spectrum
    % factor is always 2 if subsampling
    if argin.Subsample > 0
        factor = [2,0];
    else
        factor = [argin.Factor,ns-1];
    end
    % Create filter
    dFiltFFT = chFilterSpectrum(size(I),filterType,filterParams,passband,argin.Squared,factor);
    if db; imagesc(dFiltFFT); title('Detail filter'); colorbar; pause; end;
    % Apply filter
    dImageFFT = Ifft .* dFiltFFT;
    
    % *********************************
    % Project onto circular harmonics
    % *********************************
    % 0th order
    f{ns}(:,:,N+1) = argin.Weights(N+1) * real(ifft2(dImageFFT));
    if db; imagesc(f{ns}(:,:,N+1)); title('Isometric detail response'); colorbar; pause; end;
    % Other orders
    RTN = RT;
    for j = 1:N
        % Complex form
        f{ns}(:,:,N+1+j) = argin.Weights(N+1+j)*ifft2(RTN.*dImageFFT);
        if mod(j,1) == 0
            f{ns}(:,:,N+1-j) = -conj(f{ns}(:,:,N+1+j));
        else
            f{ns}(:,:,N+1-j) = conj(f{ns}(:,:,N+1+j));
        end
        f{ns}(:,:,N+1-j) = argin.Weights(N+1-j)*ifft2(conj(RTN).*dImageFFT);
        RTN = RTN .* RT;
    end
    
    % *********************************
    % Make the approximation / low-pass filter
    % *********************************
    % Approximation filter applied if subsampling or last scale
    if argin.Subsample > 0 || ns == numScales
        % Approximation filter
        factor(2) = factor(2) + 1;
        aFiltFFT = chFilterSpectrum(size(I),filterType,filterParams,'low',argin.Squared,factor);
        aImageFFT = Ifft .* aFiltFFT;
        if db; imagesc(aFiltFFT); title('Approx filter'); colorbar; pause; end;
        
        % If subsampling then downsample
        if argin.Subsample > 0
            % Downsample by 2
            I = real(ifft2(aImageFFT));
            I = I(1:2:end,1:2:end);
            Ifft = fft2(I);
        end
        if db; imagesc(I); title('Aprrox. response'); colorbar; pause; end;
        
        % If last order just save it
        if ns == numScales
            f{numScales+1}(:,:,1) = real(ifft2(aImageFFT));
        end
    end
end

end

