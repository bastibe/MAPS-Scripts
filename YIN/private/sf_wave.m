function y = sf_wave(i, samples, chans)
% y=sf_wave(i,samples,chans) - read data from file
%
%  y: data read from file (column vector or matrix)
% 
%  i: structure containing information about file
%  samples: range of samples to read ([start stop]) - default ([]) means entire file
%  chans: range of channels to read ([lo hi]) - default ([]) means all channels\
%

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

if ~nargin | ~isfield(i, 'fname')
	error('usage: i.fname = "file name", then call y=sf_wave(i)');
end
if ~isfield(i, 'format')
	i = sf_info(i);
end
% defaults
if nargin<2 | isempty(samples);
	if (i.nsamples) samples=[1 i.nsamples]; end;
end
if nargin<3 | isempty(chans)
	if (i.nchans) chans=[1 i.nchans]; end;
end

% clip
if samples(1) < 1
	samples(1) = 1;
end
if samples(2) > i.nsamples
	samples(2) = i.nsamples;
end
if samples(2) < samples(1)
	y=[];return
	error(['start sample after stop sample: ' num2str(samples(1)) '-' num2str(samples(2))]);
end
if chans(1) < 1 | chans(2) > i.nchans
	error(['requested inexistent channels: ' num2str(chans(1)) '-' num2str(chans(2))]);
end
    

% workspace matrix
if strcmp(i.format, 'matrix')
	y = i.fname(samples(1):samples(2), chans(1):chans(2));
	return
end

% use matlab functions for AU and WAV and MACSND
if strcmp(i.format, 'AU')
%	y=auread(i.fname, [samples(1) samples(2)]);
	y=audioread(i.fname, [samples(1) samples(2)]);
    y=y(:,chans(1):chans(2));
	return;
end
if strcmp(i.format, 'WAV')
%	y=wavread(i.fname, [samples(1) samples(2)]); 
	y=audioread(i.fname, [samples(1) samples(2)]); 
    y=y(:,chans(1):chans(2));
	return;
end
if strcmp(i.format, 'MACSND')
	if 3==exist('readsnd')		
		y = eval('readsnd(i.fname)');
	else
		error('cannot read MACSND on this platform');
	end
	y = y(samples(1):samples(2),chans(1):chans(2)); 
	return;
end

% close if open
% if fopen(i.fd)
% 	fclose(i.fd);
% end

if ~isfield(i, 'bytes_to_data')
	i.bytes_to_data=0;
end

% ascii formats
if strcmp(i.format, 'ascii') | strcmp(i.format, 'csv') | strcmp(i.format, 'IWAVE')
	i.fd = fopen(i.fname, 'rt');
	fseek(i.fd, i.bytes_to_data, -1);

	switch i.format
	case 'ascii'
		nsamples = samples(2) - samples(1) + 1;
		nchans = chans(2) - chans(1) + 1;
		y = zeros(nsamples, nchans);
		% skip to start
		for j=1:samples(1)-1
			line = fgetl(i.fd);
			if isempty(line) | line == -1; error('unexpected eof'); end
		end
		k=1;
		% read 
		for j=samples(1) : samples(2)
			line = fgetl(i.fd);
			if isempty(line) | line == -1; error('unexpected eof'); end
			a = sscanf(line, '%f');
			y(k,:) = a(chans(1):chans(2))';	
			k = k+1;
		end
	case 'cvs'
		error('not implemented');
	case 'IWAVE'
		error('not implemented');
	end
	fclose(i.fd);
	return
end

% binary formats
fr = samples(2) - samples(1) + 1;
skip_samples = i.nchans * (samples(1) - 1);
switch i.format
case 'uchar'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	fseek(i.fd, skip_samples * 1, 0);
	y = fread(i.fd, [fr, i.nchans], 'uchar');	
	fclose(i.fd);
case 'short'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	fseek(i.fd, skip_samples * 2, 0);
	y = fread(i.fd, [fr, i.nchans], 'short');	
	fclose(i.fd);
case 'long'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	fseek(i.fd, skip_samples * 4, 0);
	y = fread(i.fd, [fr, i.nchans], 'long');	
	fclose(i.fd);
case 'float'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	fseek(i.fd, skip_samples * 4, 0);
	y = fread(i.fd, [fr, i.nchans], 'float');	
	fclose(i.fd);
case 'double'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	fseek(i.fd, skip_samples * 8, 0);
	y = fread(i.fd, [fr, i.nchans], 'double');	
	fclose(i.fd);
case 'NIST'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	y = zeros(i.nsamples, i.nchans);
	switch i.sample_coding
	case 'pcm'
		fseek(i.fd, skip_samples * 2, 0);
		y = fread(i.fd, [fr, i.nchans], 'short');
	otherwise
		error(['cannot handle NIST sample_coding = ', i.sample_coding]);
	end
	fclose(i.fd);
case 'ESPS'
	i.fd = fopen(i.fname, 'r');
	fseek(i.fd, i.bytes_to_data, -1);
	fseek(i.fd, skip_samples * 2, 0);
	y = fread(i.fd, [fr, i.nchans], 'short');

case 'AIFF'
	i.fd = fopen(i.fname, 'r', 's');
	fseek(i.fd, i.bytes_to_data, -1);
	% should check sample size
	fseek(i.fd, skip_samples * 2, 0);
	y = fread(i.fd, [fr, i.nchans], 'short');
	fclose(i.fd);
	
otherwise
	error(['don''t know how to load format = ', i.format]);
end

y = y(:,chans(1):chans(2));
