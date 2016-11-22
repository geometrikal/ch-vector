function [ val ] = chPolyValue( coeffs, theta )
%CHPOLYVALUE - Calculate trigonometric polynomial value at a point

[r,c,N2] = size(coeffs);
N = (N2-1)/2;

kv = repmat(reshape(-N:1:N,[1,1,N2]),[r,c,1]);

rM = exp(1i * bsxfun(@times, -theta, kv));

val = real(sum(coeffs .* rM,3));

end

