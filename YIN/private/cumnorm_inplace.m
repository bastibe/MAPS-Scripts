% cumnorm_inplace(x) - cumulative mean-normalization 
%
%  x: matrix to normalize
%
% Each column of x is normalized by dividing each sample by the mean of
% samples with indices smaller or equal to it.
%
% Mex function.
% Warning: the input argument is assigned to.  This is not matlab-kosher!
