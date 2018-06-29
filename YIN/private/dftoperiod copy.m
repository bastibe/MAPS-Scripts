% prd=dftoperiod(d,b,t) - estimate period from difference function
%
%  prd: row matrix of period estimates
%
%  d: column vector or matrix of difference functions
%  b: bounds matrix ([lo, hi])
%  t: threshold
%
% Difference functions are supposed to be cumulative-mean-normalized.
% For each column of d, search for the first minimum between lo and hi that
% falls below threshold.  The index of this minimum (re 0) 
% is the period estimate.
%
% Mex function.
