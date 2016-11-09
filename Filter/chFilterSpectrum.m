function [ filterFFT, ux, uy, r, th ] = chFilterSpectrum( spectrumSize, filterType, filterParams, passband, squared, factor )
%CHFILTERSPECTRUM - Constructs the spectrum of different filters.
%
% Inputs:
%
% spectrumSize      2D size of spectrum, e.g. [512,512]
%
% filterType        Type of filter
% - lognormal:      log-normal filter
% - lognormalex:    log-normal filter with adjustable width
% - simoncelli:     Simoncelli's wavelets
% - simoncelliadj:  Simoncelli's wavelets 
% - meyer:          Meyer's first wavelet 
% - meyer2:         Meyer's second wavelet 
% - box:            Box filter
% - vow:            Variance optimised (in 1D) wavelet
% - vow2:           Variance optimised (in 2D) wavelet
% - cauchy:         Cauchy filter (like poisson)
% - papadakis:      Papadakis' wavelet
%
% filterParams      Vector of filter parameters
% - lognormal:      [wavelength,bandwidth] 
%                   bandwidth of 0.5 is wide, 0.65 medium, 0.75 narrow
% - lognormalex:    [wavelength,bandwidth,decay]
%                   decay = 2 is lognormal, higher is faster decay
% - simoncelli:     [wavelength]
% - simoncelliadj:  [wavelength]
% - meyer:          [wavelength]
% - meyer2:         [wavelength]
% - box:            [lower cutoff wavelength, higher cutoff wavelength]
% - vow:            [wavelength]
% - vow2:           [wavelength]
% - cauchy:         [wavelength,decay]
% - papadakis:      [wavelength]
%
% passband          band (default), high, low
% squared           0 - nothing, 1 - square the frequency response
% 
% factor            [factor,power] - default is [2,0]
%                   Wavelength is multiplied by factor^power. This is
%                   useful when creating wavelet sets. For example, the 1st
%                   to 3rd scales are given by [2,0], [2,1], [2,2].
%
%                   For the simoncelliadj filter, factor values other than
%                   2 will still give a partition of the spectrum. For
%                   example, one can use [sqrt(2),?] to have twice as many
%                   filters per unit log-freq unit.
%
% Outputs:
% 
% filterFFT         Spectrum of the filter
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
% Update all filters to work like simoncelliadj

% Argument parsing
if nargin < 4
    passband = 'band';
    squared = 0;
    factor = [2 0];
elseif nargin < 5
    squared = 0;
    factor = [2 0];
elseif nargin < 6
    factor = [2 0];
end

% Create FFT mesh
[ux,uy,r,th] = chFFTMesh(spectrumSize);

% Apply filter
switch lower(filterType)
    
    case {'lognormal','log-normal','loggabor','log-gabor'}
        w = filterParams(1)*factor(1)^factor(2);
        sigma = filterParams(2);
        filterFFT = exp((-(log(r*w)).^2) / (2 * log(sigma)^2));
        filterFFT(1,1) = 0;
        if strcmp(passband,'high')
            filterFFT(r >= w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r <= w) = 1;
        end
        
    case {'lognormalex','log-normalex','loggaborex','log-gaborex'}
        w = filterParams(1)*factor(1)^factor(2);
        sigma = filterParams(2);
        filterFFT = exp((-abs(log(r*w).^filterParams(3) / ...
            (filterParams(3) * (log(sigma)^filterParams(3))))));
        filterFFT(1,1) = 0;
        if strcmp(passband,'high')
            filterFFT(r >= w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r <= w) = 1;
        end
        
    case {'box','rectangle'}
        w1 = 1.0/(filterParams(1)*factor(1)^factor(2));
        w2 = 1.0/(filterParams(2)*factor(1)^factor(2));
        filterFFT = (r < w1) .* (r >= w2);
        if strcmp(passband,'high')
            filterFFT(r >= w1) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r < w2) = 1;
        end
        
    case {'cosinelog','simoncelli'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));
        filterFFT = cos(pi/2*log2(r/w));
        filterFFT(r < w/2) = 0;
        filterFFT(r > 2*w) = 0;
        if strcmp(passband,'high')
            filterFFT(r > w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r < w) = 1;
        end
    
    case {'cosinelogadj','simoncelliadj'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));
        filterFFT = cos(pi/2*log(r/w)./log(factor(1)));
        filterFFT(r < w/factor(1)) = 0;
        filterFFT(r > factor(1)*w) = 0;
        if strcmp(passband,'high')
            filterFFT(r > w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r < w) = 1;
        end
        
    case {'meyer'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));
        v1 = (2*r/w - 1); v1(v1 > 1) = 1; v1(v1 < 0) = 0;
        v2 = (r/w - 1); v2(v1 > 1) = 1; v2(v1 < 0) = 0;
        
        filterFFT = zeros(size(r));
        filterFFT(r > w/2 & r <= w) = sin(pi/2*v1(r > w/2 & r <= w));
        filterFFT(r > w & r <= 2*w) = cos(pi/2*v2(r > w & r <= 2*w));
        filterFFT(r < w/2) = 0;
        filterFFT(r > 2*w) = 0;
        if strcmp(passband,'high')
            filterFFT(r > w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r < w) = 1;
        end
        
    case {'meyer2'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));
        ra = r*(1/4)/w*8;        
        filterFFT = zeros(size(r));
        
        filterFFT(ra < 1) = 0;        
        temp = 2*sqrt(2)*(log(ra)./log(2)).^2;
        filterFFT(ra >= 1 & ra < sqrt(2)) = temp(ra >= 1 & ra < sqrt(2));        
        temp = sqrt(1 - 8*(1- log(ra)./log(2)).^4);
        filterFFT(ra >= sqrt(2) & ra < 2*sqrt(2)) = temp(ra >= sqrt(2) & ra < 2*sqrt(2));        
        temp = 2*sqrt(2)*(2 -log(ra)./log(2)).^2;
        filterFFT(ra >= 2*sqrt(2) & ra < 4) = temp(ra >= 2*sqrt(2) & ra < 4);        
        filterFFT(ra >= 4) = 0;
        
        if strcmp(passband,'high')
            filterFFT(ra >= 2) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(ra < 2) = 1;
        end
        
    case {'vow'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));        
        filterFFT = zeros(size(r));
        v1 = ...
            sqrt(1/2 + (tan(0.75*(1+2*log2(r/w))))/(2*tan(0.75)));
        v2 = ...
            sqrt(1/2 - (tan(0.75*(1+2*log2(r/(2*w)))))/(2*tan(0.75)));
        filterFFT(r > w/2 & r <= w) = v1(r > w/2 & r <= w);
        filterFFT(r > w & r <= 2*w) = v2(r > w & r <= 2*w);
        filterFFT(r < w/2) = 0;
        filterFFT(r > 2*w) = 0;
        if strcmp(passband,'high')
            filterFFT(r > w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r < w) = 1;
        end
        
    case{'vow2'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));
        
        temp1 = sqrt(log(2*r./w) ./ log(1/(4*w)));
        temp2 = sqrt(1 - log(r./w) ./ log(1/(4*w)));
        
        filterFFT = zeros(size(r));
        filterFFT(r < w/2) = 0;
        filterFFT(r >= w/2 & r < w) = temp1(r >= w/2 & r < w);
        filterFFT(r >= w & r < 2*w) = temp2(r >= w & r < 2*w);
        filterFFT(r > 2*w) = 0;
                        
        if strcmp(passband,'high')
            filterFFT(r > w) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(r < w) = 1;
        end
        
    case {'cauchy'}
        w0 = 1.0 / (filterParams(1) * factor(1)^factor(2));
        a = filterParams(2);
        sigma = a / w0;        
        nc = 1 ./ (w0.^a .* exp(- sigma * w0))^1;
        filterFFT = nc .* (r.^a .* exp(- sigma * r)).^1;
        
    case {'papa','papadakis'}
        w = 1.0/(filterParams(1)*factor(1)^factor(2));
        ra = r*(1/4)/w;        
        filterFFT = zeros(size(r));
        
        filterFFT(ra < 3/20) = 0;
        
        temp = sqrt((1 - cos(10*pi*ra - 3*pi/2))/2);
        filterFFT(ra >= 3/20 & ra < 1/4) = temp(ra >= 3/20 & ra < 1/4);
        
        filterFFT(ra >= 1/4 & ra < 3/10) = 1;
        
        temp = sqrt((1 + cos(5*pi*ra - 3*pi/2))/2);
        filterFFT(ra >= 3/10 & ra < 1/2) = temp(ra >= 3/10 & ra < 1/2);
        
        filterFFT(ra >= 1/2) = 0;
        
        if strcmp(passband,'high')
            filterFFT(ra >= 1/4) = 1;
        end
        if strcmp(passband,'low')
            filterFFT(ra < 1/4) = 1;
        end 
        
    otherwise
        error('Unknown filter type');
end

if squared == 1
    filterFFT = filterFFT.^2;
end

end

