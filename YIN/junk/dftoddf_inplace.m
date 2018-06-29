% dftoddf_inplace(d,dd,di) - double difference function from difference function
%
%  d: matrix of difference functions (time X lag)
%  dd: double difference function (lag X lag)
%  di: index (time)
%
% The double difference function at time di is calculated as a combination of
% terms of the running difference function in d. 
%
% Lags are represented implicitly by index.  First df sample is for lag=0.
%
% Mex function.
% Warning: the arguments are assigned to.  This is not matlab-kosher!

