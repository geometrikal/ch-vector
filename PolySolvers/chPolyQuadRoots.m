function [ value, theta ] = chPolyQuadRoots( beta )

% Size of 2D image RT vector matrix
[r,c,N2] = size(beta);
N = (N2-1)/2;

kv = reshape(-N:1:N,[1,1,N2]);
dbeta = bsxfun(@times,beta,kv);

ro = zeros(r,c,4);
A = dbeta(:,:,N+3);
B = dbeta(:,:,N+2);
C = dbeta(:,:,N+1);
D = dbeta(:,:,N+0);
E = dbeta(:,:,N-1);

m_min = 1e-20;
%
alpha = -3.*B.^2./(8.*A.^2) + C./A;
beta = + B.^3./(8.*A.^3) - B.*C./(2.*A.^2) + D./A;
gamma = -3.*B.^4./(256.*A.^4) + C.*B.^2./(16.*A.^3) - B.*D./(4.*A.^2) +E./A;
%
P = -alpha^2/12 - gamma;
Q = -alpha^3/108 + alpha*gamma/3 - beta^2/8;
P= -alpha.^2./12 - gamma;
Q= -alpha.^3./108 + alpha.*gamma./3 - beta.^2./8;
%
Rp= Q./2 + (Q.^2./4 + P.^3./27).^0.5;
Rm= Q./2 - (Q.^2./4 + P.^3./27).^0.5;
%
U=Rm.^(1./3);
%
mask = (U == 0);
y0=-5./6.*alpha - U;
y1=-5./6.*alpha - U + P./(3.*U);
y = y0 .* mask + y1 .* ~mask;
%
W=( alpha + 2.*y ).^(1./2);
%
ro(:,:,1) = -B./(4.*A) + 0.5 .* ( + W + ( -(3.*alpha + 2.*y + 2.*beta./(W+m_min)) ).^0.5 );
ro(:,:,2) = -B./(4.*A) + 0.5 .* ( - W + ( -(3.*alpha + 2.*y - 2.*beta./(W+m_min)) ).^0.5 );
ro(:,:,3) = -B./(4.*A) + 0.5 .* ( + W - ( -(3.*alpha + 2.*y + 2.*beta./(W+m_min)) ).^0.5 );
ro(:,:,4) = -B./(4.*A) + 0.5 .* ( - W - ( -(3.*alpha + 2.*y - 2.*beta./(W+m_min)) ).^0.5 );

theta = angle(ro);

% Evaluate at each root
val = zeros(r,c,2*N);
for j = 1:4
    val(:,:,j) = chPolyValue(beta,theta(:,:,j));
end

% Find maximum
[value, idx] = max(abs(val),[],3);

% Find theta at maximum
theta = chOpValue3rdDim(theta,idx);

end

