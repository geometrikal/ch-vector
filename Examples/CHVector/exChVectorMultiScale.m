% Load image
I = imreadgray('c_barbara.tif');

% Obtain CH vector at 4 scales
numScales = 4;
N = 3;

[f, A, theta] = chVectorMultiScale(I,N,numScales,'loggabor',[4,0.65]);

% Display results
idx = 1;
for k = 1:numScales
    for n = 0:N
        % CH vector component
        % - real part
        subplot(N+1,4,idx);
        imagesc(real(f{k}(:,:,N+1+n))); 
        daspect([1,1,1]); 
        title(['sc: ' num2str(k) ', \Re(f_' num2str(n) ')']);
        axis off;
        idx = idx+1;
        % - imag part
        subplot(N+1,4,idx);
        imagesc(imag(f{k}(:,:,N+1+n)));
        daspect([1,1,1]); 
        title(['sc: ' num2str(k) ', \Im(f_' num2str(n) ')']);
        axis off;
        idx = idx+1;
        % - amplitude
        subplot(N+1,4,idx);
        imagesc(A{k}(:,:,N+1+n));
        daspect([1,1,1]); 
        title(['sc: ' num2str(k) ', |f_' num2str(n) '|']);
        axis off;
        idx = idx+1;
        % - angle
        subplot(N+1,4,idx);
        imagesc(theta{k}(:,:,N+1+n));
        daspect([1,1,1]); 
        title(['sc: ' num2str(k) ', arg(f_' num2str(n) ')']);
        axis off;
        idx = idx+1;
    end    
    figResize(900,900);
    pause;
    idx = 1;
    close all;
end
