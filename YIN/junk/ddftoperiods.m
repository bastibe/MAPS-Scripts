% prd=ddftoperiods(d,br,bp,t) - estimate periods from double difference function
%
%  prd: 2-sample matrix of period estimates
%
%  d: column vector or matrix of difference functions
%  br: bounds matrix ([lo, hi]) for rows
%  bc: bounds matrix ([lo, hi]) for columns
%  t: threshold
%
% The double difference function is supposed to have been cumulative-mean-
% normalized along each column.
% The half-quadrant above the diagonal (row indices < column index) is searched
% for the first minimum below threshold.  The coordinates of this min give the
% period estimates.
%
% Mex function.
