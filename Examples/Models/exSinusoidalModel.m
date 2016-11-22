clear all
close all

% Load image
I = imreadgray('c_barbara.tif');

% Find the sinusoidal model at 4 scales using up to the 7th order RT
numScales = 4;
N = 7;

% Obtain the image CH vector for each scale
% - f is cell array with size numScales
% - each cell contains (rows x cols x 2N+1) matrix
f = chVectorMultiScale(I,N,numScales,'loggabor',[4,0.65]);

% Generate some weights
weights = chWeightsSinusoid(N,0,1);

% Apply to CH vector
f = chOpMultDim(f,weights);

% Obtain the model wavelet set
U = chModelVectors('sinusoidal',N,weights);

% Calculate the model (using iterative solution of the maximum)
% Create options
options = chSolveDefaultOptions
% - since we are using the sinusoidal model, the polyType is 2
options.polyType = 2;
% Solve
for k = 1:numScales
    [lambda{k}, delta{k}, orient{k}, model{k}, resid{k}] = ...
        chSolveIterative(f{k}, U, 1, options);
    
    % Amplitude
    A{k} = chOpNormDim(lambda{k},3) * sqrt(2);
    % Phase
    phase{k} = angle(lambda{k}(:,:,1) + 1i.*lambda{k}(:,:,2));
    % Resid norm
    residNorm{k} = chOpNormDim(resid{k},3);
end

% Display results
idx = 1;
for k = 1:numScales
    % Amplitude
    subplot(4,4,idx);
    imagesc(A{k}); 
    daspect([1,1,1]); 
    title(['sc: ' num2str(k) ', Amp.']);
    axis off;
    idx = idx+1;
    % Phase
    subplot(4,4,idx);
    imagesc(phase{k}); 
    daspect([1,1,1]); 
    title(['sc: ' num2str(k) ', Phase']);
    axis off;
    idx = idx+1;
    % Orientation
    subplot(4,4,idx);
    imagesc(orient{k}); 
    daspect([1,1,1]); 
    title(['sc: ' num2str(k) ', Orient.']);
    axis off;
    idx = idx+1;
     % Orientation
    subplot(4,4,idx);
    imagesc(residNorm{k}); 
    daspect([1,1,1]); 
    title(['sc: ' num2str(k) ', Resid. norm']);
    axis off;
    idx = idx+1;
end
figResize(900,900);
