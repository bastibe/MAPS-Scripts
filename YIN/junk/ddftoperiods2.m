% prd=ddftoperiods2(d,br,bp,t) - estimate periods from double difference function - specialized
%
%  prd: 2-sample matrix of period estimates
%
%  d: column vector or matrix of difference functions
%  br: bounds matrix ([lo, hi]) for rows
%  bc: bounds matrix ([lo, hi]) for columns
%  t: threshold
%
% The rectangular region between bounds is searched for the first minimum below threshold.
% The coordinates of this min give the period estimates.
%
% This version differs from ddftoperiods in that search is not limited by the diagonal.
%
% Mex function.
