function [ val, theta ] = chPolyMaxIterative( coeffs,  iterateOrders )
%CHPOLYMAXITERATIVE - Solve polynomial using iterative method Summary of this function goes here
%   Detailed explanation goes here
[r,c,N2] = size(coeffs);
N = (N2-1)/2;

iterateOrders = min(N,iterateOrders);

% Setup
od = factorial(iterateOrders);
ro = zeros(r,c,od);
pw = 2;

% Exaustive search up to iterate orders
for j = 1:iterateOrders
    % Indexes
    idx = repmat(1:j,[od/factorial(j),factorial(j)/j]);
    idx = idx(:);
    idx3 = repmat(reshape(idx,[1,1,od]),[r,c,1]);
    % Add orders
    tr = repmat(coeffs(:,:,N+1+j),[1,1,od]);
    tr = abs(tr) .* exp(1i .* (angle(tr)./j + (idx3-1)*2*pi./j));
    ro = ro + j^pw*tr;
end

% Iterate from then on
for j = iterateOrders+1:N
    % Order
    rN = coeffs(:,:,N+1+j);
    rN = abs(rN) .* exp(1i.*angle(rN)./j);
    rN = repmat(rN,[1,1,od]);
    % Angle multiplier
    am = round((angle(ro)-angle(rN)) / (2*pi/j));
    % Pick best for angle
    ro = ro + j^pw*(abs(rN) .* exp(1i.*(angle(rN)+ 2*pi*am./j)));
end
theta = mod(angle(ro),2*pi);

% Find max
kv = reshape(-N:N,1,1,N2);
for j = 1:od
    rM = exp(1i * bsxfun(@times, theta(:,:,j), kv));
    f(:,:,j) = sum(coeffs .* conj(rM),3);
end
[~, pIdx] = max(f,[],3);

% Re-organise orientation values
[ri,ci] = ndgrid(1:r,1:c);
oIdx = sub2ind([r,c,od],ri,ci,pIdx);
theta = theta(oIdx);
val = f(oIdx);

end

