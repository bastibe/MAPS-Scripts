% rdiff_inplace(x,y,d,lags,n) - in place running cross-difference function
%
%  x: column vector
%  y: column vector
%  r: result matrix (time X lag)
%  lags: 2-column matrix of lags 
%  n: (samples) frame-rate & window size (default=1)
%
% Vectors x and y are each delayed by amounts specified by lags and subtracted
% sample to sample.  The difference is squared and added over a time window
% of n samples.  This processing is repeated every n samples, as many times
% as there are columns in r.
%
% A positive lag applied to x causes it to be delayed with respect to y.
% A positive lag applied to y causes it to be delayed with respect to x.
% The number of rows of r must match that of lags. 
%
% Mex function.
% Beware: input arguments are assigned to. This is not matlab-kosher!
