function i = sf_cleanup(i)
% i=sf_cleanup(i) - cleanup after use (close i.fd if open)

% Alain de Cheveigné, CNRS/Ircam, 2002.
% Copyright (c) 2002 Centre National de la Recherche Scientifique.
%
% Permission to use, copy, modify, and distribute this software without 
% fee is hereby granted FOR RESEARCH PURPOSES only, provided that this
% copyright notice appears in all copies and in all supporting 
% documentation, and that the software is not redistributed for any 
% fee (except for a nominal shipping charge). 
%
% For any other uses of this software, in original or modified form, 
% including but not limited to consulting, production or distribution
% in whole or in part, specific prior permission must be obtained from CNRS.
% Algorithms implemented by this software may be claimed by patents owned 
% by CNRS, France Telecom, Ircam or others.
%
% The CNRS makes no representations about the suitability of this 
% software for any purpose.  It is provided "as is" without express
% or implied warranty.  Beware of the bugs.

if ~nargin ; error('sf_info: no input arguement'); end
if ~isa(i, 'struct') error('sf_info: expected struct'); end

% close file if open
if isfield(i, 'fd') & fopen(i.fd)
	fclose(i.fd);
end
