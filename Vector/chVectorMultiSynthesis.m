function [ I, Ich ] = chVectorMultiSynthesis( f, filterType, filterParams, varargin )
%CHVECTORMULTISYNTHESIS - Reconstruct image from multi-scale CH vector
%                         wavelet decomposition
%
% Inputs:
%
% f             Cells of CH vectors for each scale plus the low pass
%               component. For example f{1} to f{4} contain the image CH
%               vectors for the 1st to 4th scale and f{5} contains the low
%               pass component.
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
% I             The reconstructed image
%
% Ich           Cell array of each reconstructed scale
%
%
% Example use:  f = chVectorMultiSynthesis(f,'vow',4,'Weights',w)
%               will reconstruct the image from the CH vectors in f 
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
% Change f to a structure with all the info self-contained

% Size of image
[sr,sc,sd] = size(f{1});
% Number of CH orders
N = (sd-1)/2;
% Number of scales
numScales = numel(f)-1;

% Parse inputs
p = inputParser;
addRequired(p,'f',@iscell);
addRequired(p,'filterType');
addRequired(p,'filterParams');
addOptional(p,'Factor', 2, @isscalar);
addOptional(p,'Subsample', 0, @isscalar);
addOptional(p,'Debug', 0, @isscalar);
addOptional(p,'Squared', 0);
addOptional(p,'Weights', ones(1,2*N+1)./sqrt(2*N+1), @(x) (numel(x) == 2*N+1));

p.parse(f,numScales,filterType,filterParams,varargin{:})
argin = p.Results

db = argin.Debug;


% *********************************
% Low pass approximation
% *********************************
ns = numScales+1;
% Upsample if needed
if argin.Subsample > 0
    Ich{ns} = zeros(2*size(f{ns}));
    Ich{ns}(1:2:end,1:2:end) = f{ns};
else
    Ich{ns} = f{ns};
end
% Create approximation filter
if argin.Subsample > 0
    factor = [argin.Factor,1];
else
    factor = [argin.Factor,ns-1];
end
aFiltFFT = rtFFTFilter(size(f{ns}(:,:,1)),filterType,filterParams,'low',argin.Squared,factor);
% Apply synthesis filter
Ich{ns} = ifft2(fft2(f{ns}).*aFiltFFT);
% Add to total
I = Ich{ns};
if db; imagesc(aFiltFFT); title('Approx. filter'); colorbar; pause; end;
if db; imagesc(I); title('Approx. synthesis'); colorbar; pause; end;


% Process each scale
for ns = numScales:-1:1
    
    % Channel size
    cSize = size(f{ns}(:,:,1));
    
    % *********************************
    % Make the detail filter
    % *********************************
    % High pass if initial OR subsampling
    % Band pass if not subsampling
    if ns == 1 || argin.Subsample > 0
        passband = 'high';
    else
        passband = 'band';
    end
    % Choose factor
    % factor(1) = division of spectrum
    % factor(2) = which bandpass filter in spectrum
    % factor is always 2 if subsampling
    if argin.Subsample > 0
        factor = [2,0];
    else
        factor = [argin.Factor,ns-1];
    end
    dFiltFFT = rtFFTFilter(cSize,filterType,filterParams,passband,argin.Squared,factor);
    if db; imagesc(dFiltFFT); title('Detail filter'); colorbar; pause; end;
    
    % *********************************
    % RT filter
    % *********************************
    [ux,uy] = rtFFTMesh(cSize);
    RT = (ux + 1i*uy) ./ abs(ux + 1i*uy);
    RT(1,1) = 0;
    RTN = (ux - 1i*uy) ./ abs(ux + 1i*uy);
    RTN(1,1) = 0;
    
    % *********************************
    % Add all the orders up
    % *********************************
    % FFT of each channel
    Ffft = zeros(size(f{ns}));
    for j = 1:2*N+1
        Ffft(:,:,j) = fft2(f{ns}(:,:,j));
    end
    % Add up
    fc = zeros(size(f{ns}));
    fc(:,:,N+1) = argin.Weights(N+1).*ifft2(Ffft(:,:,N+1).*dFiltFFT);
    for j = 1:N
         n = N+1-j;
         fc(:,:,j) = argin.Weights(j).*ifft2((RT.^n).*Ffft(:,:,j).*dFiltFFT);
         fc(:,:,N+1+j) = argin.Weights(N+1+j).*ifft2((RTN.^j).*Ffft(:,:,N+1+j).*dFiltFFT);
    end
    % Add them up
    Ich{ns} = real(sum(fc,3));    
    
    % *********************************
    % Add to total
    % *********************************
    I = I + Ich{ns};
    
    % *********************************
    % Upsample if required
    % *********************************
    if ns > 1 && argin.Subsample > 0
        % Upsample
        temp = zeros(2*size(I));
        temp(1:2:end,1:2:end) = I*4;
        % Create approximation filter
        factor = [argin.Factor,1];
        aFiltFFT = rtFFTFilter(size(temp),filterType,filterParams,'low',argin.Squared,factor);
        % Filter
        I = ifft2(fft2(temp).*aFiltFFT);
    end
    
    if db; imagesc(Ich{ns}); title('Channel synthesis'); colorbar; pause; end;
    if db; imagesc(I); title('Total synthesis'); colorbar; pause; end;
end

end

