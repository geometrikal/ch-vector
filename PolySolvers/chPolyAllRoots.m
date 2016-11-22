function [ value, theta, vv ] = rtPolyAllRoots( beta )

% Size of 2D image RT vector matrix
[r,c,N2] = size(beta);
N = (N2-1)/2;

vv = zeros(r,c);

kv = reshape(-N:1:N,[1,1,N2]);
dbeta = bsxfun(@times,beta,1i*kv);
ddp = bsxfun(@times,dbeta,1i.*kv);

% Find the roots of the derivative polynomial
ro = nan(r,c,2*N);
dbeta = single(dbeta);
parfor y = 1:r
    if mod(y,25) == 0
        disp([num2str(y) '/' num2str(r) ' rows processed']);
    end
    for x = 1:c
            fd = dbeta(y,x,:);          % Pixel coefficients
            fd = fd(:);                 % Make flat
            rt = roots(fd);             % Find roots
            
            
%             fd = fd(:).';
%             n = size(fd,2);
%             r = zeros(0,1,class(fd));  
% 
%             inz = find(abs(fd) < 0.000001)
% 
%             % Strip leading zeros and throw away.  
%             % Strip trailing zeros, but remember them as roots at zero.
%             nnz = length(inz);
%             fd = fd(inz(1):inz(nnz));
%             rt = zeros(n-inz(nnz),1,class(fd));  
% 
%             % Prevent relatively small leading coefficients from introducing Inf
%             % by removing them.
%             d = fd(2:end)./fd(1);
%             while any(isinf(d))
%                 fd = fd(2:end);
%                 d = fd(2:end)./fd(1);
%             end
% 
%             % Polynomial roots via a companion matrix
%             n = length(fd);
%             if n > 1
%                 a = diag(ones(1,n-2,class(fd)),-1);
%                 a(1,:) = -d;
%                 rt = [rt;eig(a)];
%             end
%             
            
            
            
            
            tmp = nan(1,2*N);         % Store in uniform size variable
            vv(y,x) = numel(rt);
            tmp(1:numel(rt)) = rt;
            ro(y,x,:) = tmp;           % Store result
%             
%                         fd = dbeta(y,x,:);          % Pixel coefficients
%             fd = fd(:);                 % Make flat            
%             d = fd(2:end)./fd(1);
%             while any(d > 1e+5)
%                 fd = fd(2:end-1);
%                 d = fd(2:end-1)./fd(1);
%             end
%             % Polynomial roots via a companion matrix
%             n = length(d)+1;
%             if n > 1
%                 a = diag(ones(1,n-2,class(fd)),-1);
%                 a(1,:) = -d;
%                 try
%                 rt = eig(a);
%                 catch
%                     beta(y,x,:)
%                     dbeta(y,x,:)
%                     a
%                     error('adfasdf')
%                 end
%             end
%             tmp = zeros(1,2*N);         % Store in uniform size variable
%             tmp(1:numel(rt)) = rt;
%             ro(y,x,:) = tmp;            % Store result
    end
end



theta = angle(ro);
ro(abs(abs(ro)-1) > 0.01) = nan;

% Evaluate at each root
val = zeros(r,c,2*N);
for j = 1:2*N
    temp = rtPolyValue(beta,theta(:,:,j));
    %temp(abs(rtPolyValue(dbeta,theta(:,:,j))) > 0.1) = 0;
    %temp(rtPolyValue(ddp,theta(:,:,j)) > 0) = 0;
    val(:,:,j) = temp;
    %val(:,:,j) = rtPolyValue(beta,theta(:,:,j));
end
% 
% ppp = ro(136,141,:); ppp(:)
% 
% value = val;

val(isnan(ro)) = 0;

%Find maximum
[value, idx] = sort(val,3,'descend');


% Find theta at maximum
[ri,ci] = ndgrid(1:size(theta,1),1:size(theta,2),1:size(theta,3));
theta = theta(sub2ind(size(theta),ri,ci,idx));


end

