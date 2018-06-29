function i = sf_format(i)
% i=sf_format(i) - try to guess file format from 
% magic words or file suffix

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

if ~nargin; error('no argument'); end
if ~isa(i, 'struct')
	i.fname=i;
	i=sf_format(i);
	return
end
if ~isfield(i, 'fname'); error('no fname field'); end

if isa(i.fname, 'double')
	i.format = 'matrix'; return
end

if isa(i.fname, 'char')  
	disp(['guessing format of ', i.fname, '...']);
else 
	disp('guessing format...');
end

% Well-known file extensions
a = findstr('.', i.fname);
if a
	suff = i.fname( a(size(a,2))+1 : size(i.fname,2));
	suff = lower(suff);
	switch suff
	case {'short', 'long', 'float', 'double', 'ascii', 'mat'}
		i.format = suff;
		return
% the following are commented out to exercise header magic
% 	case {'aiff'}
% 		i.format = 'AIFF';
% 		return
% 	case {'aifc'}
% 		i.format = 'AIFC';
% 		return
% 	case 'wav'
% 		i.format = 'WAV';
% 		return
% 	case 'au'
% 		i.format = 'AU';
% 	case 'sdif'
% 		i.format = 'sdif';
	% others ?
	end
end

% close file if open
if exist('i.fd') & fopen(i.fd)
	fclose(i.fd);
end


% open little-endian, check for .wav magic
[i.fd, msg] = fopen(i.fname, 'r', 'l');
if i.fd == -1
	error(['could not open: >', i.fname, '<']);
	if ~isempty(msg)
		error(msg);
	end
end

% file size in bytes
if (-1 == fseek(i.fd, 0, 1)) ; error ('fseek failed');  end;	
i.nbytes = ftell(i.fd);
if i.nbytes == 0 & strcmp(computer, 'MAC2')
    disp('file is empty: perhaps its a macintosh SND resource?');
	i.format = 'MACSND';	% macintosh sound resource (maybe...)
	return
end
if i.nbytes < 4; error('file is less than 4 bytes long'); end;	
fseek(i.fd, 0, -1);			
magic = char(fread(i.fd, 4, 'uchar'))';
if strcmp(magic, 'RIFF') & i.nbytes >=12
	dummy = fread(i.fd, 1, 'ulong');	% chunk size
	magic = char(fread(i.fd, 4, 'uchar'))';
	if strcmp(magic, 'WAVE');
		i.format='WAV';
		return
	end
end

% reopen for standard binary, look for magic
fclose(i.fd);
i.fd = fopen(i.fname, 'r');
if (-1 == fseek(i.fd, 0, 1)) ; error ('fseek failed');  end;	
i.nbytes = ftell(i.fd);
fseek(i.fd, 0, -1);	
% magic strings?
magic = char(fread(i.fd, 4, 'uchar'))';
switch magic
case '.snd'
	i.format = 'AU';
	return;	
case 'FORM'
	fseek(i.fd, 8, -1);
	i.format = char(fread(i.fd, 4, 'uchar'))';
	if ~strcmp(i.format, 'AIFF') & ~strcmp(char(i.format),'AIFC')
		error (['expected AIFF or AIFC, found ' i.format]); 
	end;
	return;
case 'NIST'
	fseek(i.fd, 0, -1);
	line = fscanf(i.fd, '%s' , 1);
	if ~strcmp(line, 'NIST_1A'); 
		error(['expected NIST_1A, found ', magic]); 
	end;
	disp(i.format);	
	return;
case 'SDIF'
	i.format = 'SDIF';
	return;
end
% wff magic number?
wff_magic =  195894762;
fseek(i.fd, 0, -1);
magic = fread(i.fd, 1, 'long');
if magic == wff_magic
	i.format = 'wff';	% Nottingham IHR wff format
	return
end
% ESPS?
esps_min_size = 333; esps_magic = 27162;
if i.nbytes > esps_min_size	
	fseek(i.fd, 16, -1);
	magic = fread(i.fd, 1, 'uint');
	if magic == esps_magic
		i.format = 'ESPS';
		return
	end
end


% Magic didn't succeed, try to determine type based on file name
a = findstr('.', i.fname);
if a
	suff = i.fname( a(size(a,2))+1 : size(i.fname,2));
	switch suff
	case 's'
		i.format = 'short';
		return
	case 'f'
		i.format = 'float';
		return
		% others ?
	end
end


% Ascii formats
fclose(i.fd);
i.fd = fopen(i.fname, 'rt');
line = fscanf(i.fd, '%s' , 1);
if strcmp(line, '|WAVE'); 
	i.format = 'iwave';			% Nottingham IHR |WAVE format
	return;
end

% plain ascii ?
limit = 2048;	% test only this many bytes
hunk = char(fread(i.fd, min(limit, i.nbytes)));
chars = unique(hunk);
if isempty(setdiff(chars, [9, 10, 13, ...
	double(' abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890+-_.,:;')]))
	% ascii
	if ~isempty(setdiff(chars, [9, 10, 13, double(' 1234567890eE+-.,')]))
		% not just numerical: try removing first line
		fseek(i.fd, 0, -1);
		line = fgets(i.fd);
		if -1 == line; return ; end	% line too long 
		i.bytes_to_data = ftell(i.fd);
		hunk = char(fread(i.fd, min(limit, i.nbytes-ftell(i.fd))));
		chars = unique(hunk);
		warning('skipping first line of what seems like an ascii file');
	end 
	if isempty(setdiff(chars, [9, 10, 13, double(' 1234567890eE+-.')]))
		i.format = 'ascii';
		return
	end
	if isempty(setdiff(chars, [9, 10, 13, double(' 1234567890eE+-.,')]))
		w.format = 'csv';
		return
	end
end

% Default
if i.nbytes == 2 * round(i.nbytes/2);	% even number of bytes
	i.format = 'short';
else
	i.format = 'uchar';
end
			
