function [ lambda, delta, orient, model, resid ] = chSolveIterative( f, U, K, options )
%CHSOLVEITERATIVE - Finds minumum distance for a wavelet set, then repeats 
%K times with the residual
%
% Inputs:
%
% f             Weighted image CH vector matrix (rows x cols x 2N+1)
%               (N is maximum RT order)
% U             Weighted model matrix of CH vectors (M x 2N+1)
%               (M is number of vectors of model)
% K             Number of iterations to perform
%
% options       Options is a structure containing the following fields:
%   method        Method of calculating distance polynomial maximum
%   - roots       Use roots of derivative polynomial
%                 params: no parameters
%   - brute       Evaluate at fixed points, choose best point, refine
%                 params: [number of points,number of refinements]
%   * iterative   Estimate minimum by permutation estimate then iteration
%                 params: [minimum number of orders before iteration]
%                 A value of 2* or 3 is recommended
%
%   params        Parameters of method as above
%   polyType      Effective order of model polynomial. For example, the
%                 half-polyTypeal model is 1, but the sinusoidal model is 2
%                 as orientation is only defined over 2\pi/2 (180 deg)
%   regType       Use regularisation when solving
%   * none        Default
%   - l2          l2 regularisation
%   - l1approx    Psuedo-l1 that generates solutions using no reg. then
%                 applies an l1 penalty
%   regPenalty    Regularisation penalty amount
%   positive      0* - nothing, 1 - restrict model amplitudes to positive
%                 values. Only valid for single wavelet models.
%
% All options fields are optional. If left out, defaults will be used.
%
%
% Outputs:    
%
% lambda        Amplitude of wavelet in model (rows x cols x M x K)
% delta         Amplitude of psuedo-inverse wavelet in model (as above)
% orient        Orientation for each iteration (rows x cols x K)
% model         Model vector for each iteration (rows x cols x 2N+1 x K)
% resid         Model vector for each iteration (rows x cols x 2N+1 x K)
%
%
% Written by:
%
% Ross Marchant
% James Cook University
% ross.marchant@my.jcu.edu.au
%

% Default options
if nargin < 4
    options.method = 'iterative';
    options.params = 0;
end

if ~isfield(options, 'polyType')
    options.polyType = 0;
end

if ~isfield(options, 'smoothing')
    options.smoothing = 0;
end

if ~isfield(options, 'positive')
    options.positive = 0;
end

if ~isfield(options, 'regType')
    options.regType = 'none';
end

if ~isfield(options, 'regPenalty')
    options.regType = 'none';
end


% Check inputs
if ndims(U) > 2
    error('Wavelet CH matrix must only have 1 dimension');
end
if max(size(U)) ~= size(f,3)
    error('Wavelet CH matrix must be M x 2N+1 where M is number of wavelets');
end

% Transpose U if necessary
if size(U,1) < size(U,2)
    U = U.';
end

% Regulisation penalty
alpha = options.regPenalty;


% Get size of image (r,c), max RT order (N), and wavelets in model (M)
[r,c,N2] = size(f);
N = (N2-1)/2;
M = size(U,2);

% Calculate psuedo-inverse
if ~strcmpi(options.regType,'l2')
    % - normal
    Uplus = pinv(U);
else
    % - with l2 penalty
    Uplus = (U'*U + alpha*eye(M))\U';    
end

% Pre-allocate
model = zeros([size(f),K]);
resid = zeros([size(f),K]);


disp('Solving model...')
if ~options.positive || M > 1
    disp('- amplitudes will be both positive and negative');
else
    disp('- amplitudes will only be positive');
end

% Loop through each iteration
for k = 1:K
    disp(['* iteration: ' num2str(k) ' of ' num2str(K)]);
    
    % Solve for positive and negative amplitudes?    
    if ~options.positive || M > 1
        
        % Start timer
        tstart = tic;
        
        % Pre-allocate memory
        pch = zeros(size(f,1),size(f,2),4*N+1,M);
        deltaP = zeros(size(f,1),size(f,2),2*N+1,M);
        lambdaP = zeros(size(f,1),size(f,2),2*N+1,M);
        
        % Calculate reponse polynomial
        for j = 1:M
            deltaP(:,:,:,j) = bsxfun(@times,f,permute(conj(U(:,j)),[3,2,1]));
            lambdaP(:,:,:,j) = bsxfun(@times,f,permute(Uplus(j,:),[3,1,2]));
            pch(:,:,:,j) = chOpConv3rdDim(lambdaP(:,:,:,j),deltaP(:,:,:,j));
            if strcmpi(options.regType,'l2')
                pch(:,:,:,j) = pch(:,:,:,j) - alpha * chOpConv3rdDim(lambdaP(:,:,:,j),lambdaP(:,:,:,j));
            end
        end
        p = sum(pch,4);
        
        % Remove zero orders?
        timing(1) = toc(tstart);
        if options.polyType == 2
            p = p(:,:,1:2:end);
        elseif options.polyType == 3
            p = p(:,:,mod(-2*N:2*N,3) == 0);
        elseif options.polyType == 4
            p = p(:,:,mod(-2*N:2*N,4) == 0);
        end
        
        % Smooth polynomial co-efficients?
        if options.smoothing ~= 0
            for si = 1:size(p,3)
                p(:,:,si) = imgaussian(p(:,:,si),options.smoothing,options.smoothing*6);
            end
        end
        
        % Show polynomial at debug location?    
        if isfield(options, 'debug')
            for j = 1:M
                pp = deltaP(options.debug(1),options.debug(2),:,j);
                pp = pp(:);
                chPlotTrigPoly(pp,-1,100); title('delta'); pause;
                pp = lambdaP(options.debug(1),options.debug(2),:,j);
                pp = pp(:);
                chPlotTrigPoly(pp,-1,100); title('lambda'); pause;
                pp = pch(options.debug(1),options.debug(2),:,j);
                pp = pp(:);
                chPlotTrigPoly(pp,-1,100); title('combined'); pause;
            end
            pp = p(options.debug(1),options.debug(2),:);
            pp = pp(:);
            rtPlotTrigPoly(pp,-1,100); pause;
        end
        
    % Solve for positive amplitudes only?
    else
        lambdaP(:,:,:,1) = bsxfun(@times,f,permute(conj(U(:,1)),[3,2,1]));
        p = lambdaP(:,:,:,1);
        if options.polyType > 1
            p = p(:,:,mod(-N:N,options.polyType) == 0);
        end
    end
    
    % Find maximum using roots method?
    if strcmpi(options.method,'roots')
        
        % Apply pseudo-l1 approximation?  
        if strcmpi(options.regType,'l1approx')             
                    
            [val, or] = chPolyAllRoots(p);            
            sumv = zeros(size(val));
            for oi = 1:size(or,3)
                for chi = 1:M
                    sumv(:,:,oi) = sumv(:,:,oi) - alpha.*abs(chPolyValue(lambdaP(:,:,:,chi),or(:,:,oi)));
                end
            end
            val = val + sumv;            
            [~,idx] = max(val,[],3);
            orient(:,:,k) = rtValue3rdDim(or,idx);
            
        % Only solve for positive amplitudes?
        elseif options.positive || M == 1            
            [val, orient(:,:,k)] = chPolyMaxPositiveRoots(p, deltaP, lambdaP);
            
        % Solve special case using quadratic (fast)
        elseif options.polyType == 2 && M == 2
            [val, orient(:,:,k)] = chPolyQuadRoots(p);
        
        % Solve as normal
        else            
            tstart = tic;
            [val, orient(:,:,k)] = chPolyMaxRoots(p);
            timing(2) = toc(tstart);
        end
        
    % Solve using iterative method?
    elseif strcmpi(options.method,'iterative')
        tstart = tic;        
        [val, orient(:,:,k)] = chPolyMaxIterative(p,options.params(1));
        timing(2) = toc(tstart);
        
    % Otherwise error
    else
        error(['The ' options.method ' is not currently supported'])
    end
    
    % Get modulo of orientation value if poly has some zero orders
    if options.polyType > 1
        orient(:,:,k) = mod(orient(:,:,k)/options.polyType,2*pi/options.polyType);
    end
    
    % Find amplitude of each wavelet
    for j = 1:M
        lambda(:,:,j,k) = chPolyValue(lambdaP(:,:,:,j), orient(:,:,k));
        if ~options.positive && size(U,2) > 1
            delta(:,:,j,k) = chPolyValue(deltaP(:,:,:,j), orient(:,:,k));
        else
            delta = 0;
        end
    end
     
    % Create model
    for j = 1:M
        chan = permute(U(:,j),[3,2,1]);
        chan = repmat(chan,r,c,1);
        chan = chVectorRotate(chan,orient(:,:,k),3);
        chan = bsxfun(@times,chan,lambda(:,:,j,k));
        model(:,:,:,k) = model(:,:,:,k) + chan;
    end
    
    % Create residual
    resid(:,:,:,k) = f - model(:,:,:,k);
    f = resid(:,:,:,k);
end

end

