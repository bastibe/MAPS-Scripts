% rsum_inplace(x,N) - in place running sum
%
%  x: matrix to sum
%  N: samples - size of running sum window
%
% Matrix x is summed row-wise.  The last (N-1) samples are not valid.
% N may be fractionary (handled by linear interpolation).
% 
% Mex function.  
% Beware: input argument is assigned to.  This is not matlab-kosher!
