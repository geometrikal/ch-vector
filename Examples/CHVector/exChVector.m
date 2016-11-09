% Load image
I = imreadgray('c_barbara.tif');

% Obtain CH vector
N = 3;

[f, A, theta] = chVector(I,N,'loggabor',[32,0.65]);

% Display results
idx = 1;
for n = 0:N
    % CH vector component
    % - real part
    subplot(N+1,4,idx);
    imagesc(real(f(:,:,N+1+n))); 
    daspect([1,1,1]); 
    title(['\Re(f_' num2str(n) ')']);
    axis off;
    idx = idx+1;
    % - imag part
    subplot(N+1,4,idx);
    imagesc(imag(f(:,:,N+1+n)));
    daspect([1,1,1]); 
    title(['\Im(f_' num2str(n) ')']);
    axis off;
    idx = idx+1;
    % - amplitude
    subplot(N+1,4,idx);
    imagesc(A(:,:,N+1+n));
    daspect([1,1,1]); 
    title(['|f_' num2str(n) '|']);
    axis off;
    idx = idx+1;
    % - angle
    subplot(N+1,4,idx);
    imagesc(theta(:,:,N+1+n));
    daspect([1,1,1]); 
    title(['arg(f_' num2str(n) ')']);
    axis off;
    idx = idx+1;
end    
figResize(900,900);
