clear all
close all

% Test image
[I, Iphase] = chImageSinusoid(512,32,1,0,pi/4);
Iphase = angle(cos(Iphase) + 1i.*sin(Iphase));

% Use up to the 7th order RT
N = 7;

% Image CH vector at sinusoid wavelength
f = chVector(I,N,'loggabor',[32,0.65]);

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
[lambda, delta, orient, model, resid] = chSolveIterative(f, U, 1, options);
    
% Calculate the sinusoid amplitude and phase
% Amplitude
A = chOpNormDim(lambda,3) * sqrt(2);
% Phase
phase = angle(lambda(:,:,1) + 1i.*lambda(:,:,2));
% Resid norm
residNorm = chOpNormDim(resid,3);

% Display results
idx = 1;
% Amplitude
subplot(1,5,idx);
imagesc(A); 
daspect([1,1,1]); 
title(['Amp.']);
axis off;
colorbar;
idx = idx+1;
% Phase
subplot(1,5,idx);
imagesc(phase); 
daspect([1,1,1]); 
title(['Phase']);
axis off;
colorbar;
idx = idx+1;
% Phase actual
subplot(1,5,idx);
imagesc(Iphase); 
daspect([1,1,1]); 
title(['Correct phase']);
axis off;
colorbar;
idx = idx+1;
% Orientation
subplot(1,5,idx);
imagesc(orient); 
daspect([1,1,1]); 
title(['Orient.']);
axis off;
colorbar;
idx = idx+1;
 % Residual norm
subplot(1,5,idx);
imagesc(residNorm); 
daspect([1,1,1]); 
title(['Resid. norm']);
axis off;
idx = idx+1;
colorbar;
figResize(1400,200);
figResize(1400,200);
