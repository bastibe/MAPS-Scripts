% minparabolic2_inplace(d,xoffset,yoffset) - paraboloid interpolation near minima
%
%  d: matrix to interpolate
%  xoffset: horizontal (col) offset of interpolated minimum
%  yoffset: vertical (row) offset of interpolated minimum
%
% At each local minimum of d (sample smaller or equal to its 4 immediate
% neighbors), a paraboloid with equation 
%   z = a*x^2 +b*x + c*y^2 + d*y + e
% (no x*y term) is fit to the 5 samples. The value of the center sample is
% replaced by the minimum of the paraboloid. The position of the minimum is
% stored in xoffset and yoffset as vertical and horizontal offsets relative 
% to the center sample.
%
% Mex function.
% Beware: arguments are assigned to. This is not matlab-kosher!
