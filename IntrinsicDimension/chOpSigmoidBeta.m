function [ y ] = chOpSigmoidBeta( x, halfway, steepness, clamp )
%CHOPSIGMOIDBETA Sigmoidal projection using the Beta function
% Sigmoidal function where input in the range [0,1] will always be 
% projected onto the range [0,1] such that 0 -> 0 and 1 -> 1
%
% Inputs:
% x         Input signal in range [0,1]
% halfway   Value of x where y = 0.5 (approximately)
% steepness Steepness of the transition at the halfway point approximately
%
% Outpus:
% y         Output of projection in range [0,1]

if ~exist('clamp')
    clamp = [0,1];
end

x = (x - clamp(1)) ./ diff(clamp);
x(x<0) = 0;
x(x>1) = 1;
halfway = (halfway - clamp(1)) ./ diff(clamp)

flag = 0;
if halfway > 0.5
    x = 1-x;
    halfway = 1-halfway;
    flag = 1;
end

a = steepness
b = steepness*(1/halfway - 1);
x = x;
x(x > 1) = 1;
x(x < 0) = 0;

phiadj = betainc(x,a,b);
phiadj = phiadj;

if flag == 1
    phiadj = 1-phiadj;
end

y = phiadj;

end

