function [ options ] = chSolveDefaultOptions( )
%CHSOLVEDEFAULTOPTIONS - Create default options structure

options.method = 'iterative';
options.params = 3;
options.polyType = 0;
options.positive = 0;
options.regType = 'none';
options.regPenalty = 0;

end

