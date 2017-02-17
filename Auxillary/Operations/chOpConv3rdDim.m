function [ r ] = chOpConv3rdDim( a, b )
%CHOPCONVDIM - convolve along 3rd dimension

if numel(which('inplaceprod')) == 0
    error(['Please run convnfft_install.m from within the Auxillary/External/convnfft directory' ...
           ', or change alg = 1 to alg = 0 in chOpConvDim.m']);
end

if nargin < 2
    issame = 1;
else
    issame = 0;
end

alg = 1;

if(size(a,3) == 1)
    r = a.*b;
    disp('11')
    return;
end

if alg == 1
    if issame
        as = size(a)
        a = reshape(a,as(1)*as(2),as(3)).';
        a = padarray(a,[(as(3)-1)/2,0],0,'both');
        a = fft(a);
        r = ifft(a.*a);
        r = reshape(r,as(1),as(2),size(r,1));
    else   
        as = size(a);
        bs = size(b);
        a = reshape(a,as(1)*as(2),as(3));
        b = reshape(b,bs(1)*bs(2),bs(3));
        r = convnfft(a,b,'full',2);
        r = reshape(r,as(1),as(2),size(r,2));
    end      
else
    sa = size(a);
    sb = size(b);
    % Pad a
    aPad = padarray(a,[0,0,sb(3)-1],0,'both');
    % Flip b
    b = b(:,:,end:-1:1);
    % Result
    r = zeros(sa(1),sa(2),sa(3)+sb(3)-1);
    % Convolve
    sz = sa(3)+sb(3)-1;
    for k = 1:ceil(sz/2)
        r(:,:,k) = sum(aPad(:,:,k:sb(3)-1+k).*b,3);
    end
    c=1;
    for k = ceil(sz/2)+1:sz
        r(:,:,k) = conj(r(:,:,k-2*c));
        c=c+1;
    end
end

end

