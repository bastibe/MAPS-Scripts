% [y,idx] = minparabolic(x) - locate minima using parabolic interpolation
%
%  x: input vector
%  y: same as x, with local minima replaced by interpolated minima
%  idx: position of interpolated minima (fractionary)
%
% At each local minimum of x, that sample and its neighbors are used to
% determine a parabola.  The minimum of the parabola is used as the new
% minimum value, the abscissa of the minimum is used as its position.
