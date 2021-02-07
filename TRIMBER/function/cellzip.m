function [ziped] = cellzip(f,a,b)
% CELLZIP  Zip two cell arrays by a function
%
%   [ZIPED] = CELLZIP(F,A,B)
%
%   Computes ZIPED{i} = F(A{i},B{i}) foreach i in A,B.

ziped = arrayfun(@(i) f(a{i},b{i}),1:length(a),'Uniform',false);
