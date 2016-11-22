function [ value, theta ] = rtPolyMaxRoots( beta )

% Size of 2D image RT vector matrix
[r,c,N2] = size(beta);
N = (N2-1)/2;

kv = reshape(-N:1:N,[1,1,N2]);
dbeta = bsxfun(@times,beta,kv);

% Find the roots of the derivative polynomial
ro = zeros(r,c,2*N);
dbeta = single(dbeta);
parfor y = 1:r
    if mod(y,25) == 0
        disp([num2str(y) '/' num2str(r) ' rows processed']);
    end
    for x = 1:c
            fd = dbeta(y,x,:);          % Pixel coefficients
            fd = fd(:);                 % Make flat
            rt = roots(fd);             % Find roots
            tmp = zeros(1,2*N);         % Store in uniform size variable
            tmp(1:numel(rt)) = rt;
            ro(y,x,:) = tmp;            % Store result
    end
end
% Convert result into an angle
theta = angle(ro);

% Evaluate at each root
val = zeros(r,c,2*N);
for j = 1:2*N
    val(:,:,j) = rtPolyValue(beta,theta(:,:,j));
end

% Find maximum
[value, idx] = max(abs(val),[],3);

% Find theta at maximum
theta = rtValue3rdDim(theta,idx);

end

