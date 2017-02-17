function [ g ] = chVectorRotate( f, theta, dimension )
%CHVECTORROTATE - Rotate the CH vector

r = chRotationVector(theta, (size(f,dimension)-1)/2, dimension);
g = bsxfun(@times,f,r);

end

