function [ g ] = chRotateVector( f, theta, dimension )
%CHVECTORROTATE - Rotate the CH vector
size(f)
r = chRotationVector(theta, (size(f,dimension)-1)/2, dimension);
g = bsxfun(@times,f,r);

end

