function i = sf_info(i)
% i=sf_info(i) - extract useful info from file

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
if ~isa(i, 'struct')
	j.fname=i;
    i=j;
	i=sf_info(i);
	return
end
if ~isfield(i, 'fname'); error('sf_info: no fname field'); end

% guess format only if unknown
if ~isfield(i, 'format') || isempty(i.format);
	i = sf_format(i);
	if ~strcmp(i.format, 'matrix') disp(i.format); end
end

% handle workspace matrices as if they were files
if strcmp(i.format, 'matrix')
	[nrows, ncols] = size(i.fname);
	i.nchans=ncols;
	i.nsamples=nrows;
	i.totalsamples = i.nsamples*i.nchans;
    i.sr=[];
	return;
end

% close file if open
if exist('i.fd') & fopen(i.fd)
	fclose(i.fd);
end

% use standard matlab functions for AU and WAV and MACSND
if strcmp(i.format, 'AU')
	if isempty(findstr('.', i.fname))
		disp(['WARNING: matlab function AUREAD requires '...
		'.au or .snd suffix on file name: ' i.fname]);
	end
%	sz = auread(i.fname, 'size');
    [x,sr]=audioread(i.fname); sz=size(x);
	i.nsamples=sz(1);
	i.nchans=sz(2);
	i.totalsamples = i.nsamples*i.nchans;
    i.sr=sr;
	%[dummy, i.sr, i.samplebits] = auread(i.fname, 1);
	return;
end
if strcmp(i.format, 'WAV')
	if isempty(findstr('.wav', i.fname)) & isempty(findstr('.WAV', i.fname))
		disp(['WARNING: expected '...
		'.wav suffix on file name: ' i.fname]);
	end
	%sz = wavread(i.fname, 'size');
    [x,sr]=audioread(i.fname); sz=size(x);
	i.nsamples=sz(1);
	i.nchans=sz(2);
	i.totalsamples = i.nsamples*i.nchans;
	%[dummy, i.sr, i.samplebits] = wavread(i.fname, 1);
    i.sr=sr;
	return;
end
if strcmp(i.format, 'MACSND')
	% must load the data to get info - this is stupid
	if ~isempty(findstr(':', i.fname))
		disp(['matlab function READSND cannot handle an '...
		'indirect path: ' i.fname]);
	end
	if 3==exist('readsnd')
		[data, i.sr] = eval('readsnd(i.fname)');
	else
		error('cannot read MACSND on this platform');
	end
	i.nsamples = size(data,2);
	i.nchans = size(data,1);
	i.totalsamples = i.nsamples*i.nchans;
	return
end

% reopen file
if strcmp(i.format, 'ascii') ...
	| strcmp(i.format, 'csv') | strcmp(i.format, '|WAVE')
	[i.fd, msg] = fopen(i.fname, 'rt');
else
	[i.fd, msg] = fopen(i.fname, 'r', 'ieee-be.l64');
end
if i.fd == -1
	if isempty(msg)
		error(['could not open: ', i.fname]);
	else
		error(msg)
	end
end
fd = i.fd;
if ~isfield(i, 'nbytes')
	if (-1 == fseek(i.fd, 0, 1)) ; error ('fseek failed');  end;
	i.nbytes = ftell(i.fd);
end
fseek(fd, 0, -1);	% rewind


switch i.format
%%%%%%%%%%%%%%%%%%% AIFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'AIFF','AIFC'}
	fseek(fd, 12, 0); % skip container chunk
	% skip over spurious chunks
        idx=ftell(fd);
	while 1
        	magic=char(fread(fd,4,'uchar'))';
                if strcmp(magic,'COMM'); break; end;
                idx = idx+1;
                status = fseek(fd,idx,-1);
		if status == -1
			error('expected COMM magic word, found eof');
		end;
	end;

        %ckSize=fread(fd,1,'int32');
	%status = fseek(fd, ckSize, 0);		% skip to end of chunk
        %if status == -1
        %      error('unexpected eof');
        %end

	%while (1)
	%	magic = char(fread(fd, 4, 'char'))';
	%	if ~strcmp(magic, 'SSND')
	%		fseek(fd, -4, 0);	% skip back
	%		break;
	%	end
	%	ckSize = fread(fd, 1, 'long');
	%	fseek(fd, ckSize, 0);	% skip to end of sound chunk
	%end
	%magic = char(fread (fd, 4, 'char'))';
	%if ~strcmp(magic, 'COMM'); error(['expected COMM, found ', magic]) ; end
	commsz = fread(fd, 1, 'int32');
	i.nchans = fread(fd, 1, 'int16');
	i.nsamples = fread(fd, 1, 'uint32');
	i.totalsamples = i.nsamples*i.nchans;
	i.samplebits = fread(fd, 1, 'int16');
	switch i.samplebits
	case 16
		i.sample_bytes = 2;
		i.sample_type = 'int16';
	case 32
		i.sample_bytes = 4;
		i.sample_type = 'int32';	% or float?
	otherwise
		error(['unexpected samplebits: ' num2str(i.samplebits) ]);
	end
	% read sampling rate using Hideki Kawahara's code:
	srex1=fread(fd,1,'uint16');
	srex2=fread(fd,1,'uint64');
	if strcmp(char(i.format),'AIFC')
	    compress=fread(fd,4,'uchar');
		if ~strcmp(char(compress),'NONE')
		    error('Compression is not supported.');
		end;
	    	fseek(fd, commsz-22, 0);
	end;
	i.sr = 2^(srex1-16383)*srex2/hex2dec('8000000000000000');
	%fseek(fd, 12, -1);	% skip back to end of container chunk
	% skip over eventual common chunk
	%while(1)
	%	magic = char(fread(fd, 4, 'char'))';
	%	if ~strcmp(magic, 'COMM')
	%		fseek(fd, -4, 0);	% skip back
	%		break;
	%	end
	%	ckSize = fread(fd, 1, 'long');
	%	fseek(fd, ckSize, 0);	% skip over chunk
	%end
	magic=char(fread(fd,4,'uchar'))';
	while ~strcmp(char(magic),'SSND')
		[ckSize, count]=fread(fd,1,'int32');
		if ~count;
			error('expected chunk size field, found eof');
			return
		end
		status = fseek(fd, ckSize, 0);		% skip to end of chunk
		if status == -1
			error('expected SSND magic word, found eof');
			return
		end;
		magic=char(fread(fd,4,'uchar'))';
	end;

	%magic = char(fread(fd, 4, 'char'))';
	%if ~strcmp(magic, 'SSND')
	%	error (['expected SSND, found' magic]);
	%end
	fseek(fd, 12, 0);	% skip over ckSize, offset and blocksize fields
	i.bytes_to_data = ftell(fd);
	if i.totalsamples*i.sample_bytes ~= i.nbytes-i.bytes_to_data
		disp(['WARNING: header fields sample_bytes: ' ...
		num2str(i.sample_bytes)]);
		disp (['and sample and channel count: ' num2str(i.nsamples) ...
		', ' num2str(i.nchans)]);
		disp(['are inconsistent with offset to data: ' ...
		num2str(i.bytes_to_data) ' and file size: ' ...
		num2str(i.nbytes)]);
		disp (['(' num2str(i.totalsamples*i.sample_bytes) ...
		' ~= ' num2str(i.nbytes-i.bytes_to_data) ')']);
	end

	return
%%%%%%%%%%%%%%%%%%% NIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'NIST'
%	fseek(fd, 0, -1);
%	line = fscanf(fd, '%s' , 1);
%	if ~strcmp(line, 'NIST_1A'); error(['expected NIST_1A, found ', magic]); end
	fseek(fd, 8, 0);	% skip over magic string
	i.bytes_to_data = fscanf(fd, '%d', 1);
	while (1)
		key = fscanf(fd, '%s', 1);
		if strcmp(key, 'end_head'); break; end;
		% read third field according to spec in second field (type)
		% this may need refining...
		type = fscanf(fd, '%s', 1);
		if strcmp(type(1:2), '-s')
			bytes_to_read = sscanf(type(3:end), '%d', 1);
			fseek(fd, 1, 0);	% skip blank
			value = char(fread(fd, bytes_to_read, 'char'))';
		else
			value = fscanf(fd, '%f', 1);
		end
		i = setfield(i, key, value);
	end
	% give standard names to useful fields
	if isfield(i, 'channel_count'); i.nchans = i.channel_count; end
	if isfield(i, 'sample_count'); i.nsamples = i.sample_count; end
	i.totalsamples = i.nsamples*i.nchans;
	if isfield(i, 'sample_rate');
		i.sr = i.sample_rate;
		i.xunits = 's';
	end
	if ~isfield(i, 'sample_coding'); i.sample_coding = 'pcm'; end
	i.bytes_to_data = 1024; 	% needs checking
	i.sample_bytes=2;
	i.sample_type='int16';
	return
%%%%%%%%%%%%%%%%%%% |WAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case '|WAV'
	line = fscanf(fd, '%s' , 1);	% skip first line
	i.nsamples = fscanf(fd, '%d', 1);
	i.nchans = fscanf(fd, '%d', 1);
	i.totalsamples = i.nsamples*i.nchans;
	i.bytes_to_data = ftell(fd);
	% channel info handled during data read
	i.totalsamples = i.nsamples * i.nchans;
	return
%%%%%%%%%%%%%%%%%%% WFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'wff'
	fseek(fd, 4, 0);	% skip magic number
	i.version = fread(fd, 1, 'long');
	i.type_code = fread(fd, 1, 'long');
	i.nchans = fread(fd, 1, 'long');
	i.info.channel_flags = fread(fd, 1, 'long');
	i.bytes_to_data = fread(fd, 1, 'long');
	fseek(fd, 40, 0);
	i.gen_prog_name = fread(fd, 32, 'char');
	i.comment = fread(fd, 32, 'char');
	i.sample_bytes = 2;
	i.sample_type = 'int16';
	return;
% channel info handled later during channel read
%%%%%%%%%%%%%%%%%%% ESPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ESPS'
% Based on Peter Kabal's afsp package.
% This handles at least one kind of ESPS waveform file.  Others?
	i.machine_code = fread(fd, 1, 'uint');
	i.version_check_code = fread(fd, 1, 'uint');
	i.bytes_to_data = fread(fd, 1, 'uint');
	i.record_size = fread(fd, 1, 'uint');
	fseek(fd, 20, -1);
	i.EDR_ESPS_flag = fread(fd, 1, 'uint');
	i.align_pad_size = fread(fd, 1, 'uint');
	fseek(fd, 32, -1);
	i.file_type = fread(fd, 1, 'uint16');
	fseek(fd, 40, -1);
	i.file_creation_date_time = char(fread(fd, 26, 'char'))';
	i.header_version = char(fread(fd, 8, 'char'))';
	i.program_name = char(fread(fd, 16, 'char'))';
	i.program_version = char(fread(fd, 8, 'char'))';
	i.compile_date = char(fread(fd, 26, 'char'))';
	i.tag = fread(fd, 1, 'uint');
	fseek(fd, 132, -1);
	i.ndoubles = fread(fd, 1, 'uint');
	i.nfloats = fread(fd, 1, 'uint');
	i.nlongs = fread(fd, 1, 'uint');
	i.nshorts = fread(fd, 1, 'uint');
	i.nchars = fread(fd, 1, 'uint');
	i.fixed_header_size = fread(fd, 1, 'uint');
	i.var_header_size = fread(fd, 1, 'uint');
	%w.dunno_what = fread(fd, 1, 'uint');
	fseek(fd, 160, -1);
	i.user = char(fread(fd, 8, 'char'))';
	% scan the rest of the header to find sampling rate
	a = ftell(fd);
	bytes_left_in_header = i.bytes_to_data - ftell(fd);
	hunk = char(fread(fd, bytes_left_in_header, 'uchar'))';
	b = findstr(hunk, 'record_freq');
	%fseek(fd, a+b-1, -1);
	%w.xxxx = char(fread(fd, 12, 'char'))';
	%w.count = fread(fd, 1, 'uint');
	%w.data_code = fread(fd, 1, 'ushort');
	fseek(fd, a+b-1+12+2+4, -1);
	i.sr = fread(fd, 1, 'double');
	i.nchans = i.ndoubles+i.nfloats+i.nlongs ...
		+i.nshorts+i.nchars;
	recordsize = i.ndoubles*8 + i.nfloats*4 + i.nlongs*4 ...
		+i.nshorts*2+i.nchars;
	i.nsamples = (i.nbytes - i.bytes_to_data) / recordsize;
	i.totalsamples = i.nsamples * i.nchans;
	i.sample_bytes = 2;			% bug
	i.sample_type = 'int16'; 	% bug
	return;
%%%%%%%%%%%%% PLAIN ASCII %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'ascii'
	if isfield(i, 'bytes_to_data')
		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	line = fgetl(fd);
	i.nchans = size(sscanf(line, '%f'), 1);
	nlines = 1;
	while (1)
		line = fgets(i.fd);
		if isa(line, 'double') & -1 == line; break; end
		nlines = nlines+1;
	end
	i.nsamples = nlines;
	i.totalsamples = i.nsamples*i.nchans;
    i.sr=1;
	return
%%%%%%%%%%%%% COMMA-SEPARATED VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'csv';
	if isfield(i, 'bytes_to_data')
		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	% todo
	return
%%%%%%%%%%%%% Binary  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case {'uchar'}; % should take care of signed/unsigned
	i.sample_bytes = 1;
	i.sample_type = 'uchar';
	if isfield(i, 'bytes_to_data')
		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	if ~isfield(i, 'nchans')
		i.nchans = 1;
	end
	i.nsamples = (i.nbytes-i.bytes_to_data)/i.sample_bytes;
	i.totalsamples = i.nsamples * i.nchans;
    i.sr=1;
	return
case {'short', 'int16'};
	i.sample_bytes = 2;
	i.sample_type = 'int16';
	if isfield(i, 'bytes_to_data')
		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	if ~isfield(i, 'nchans')
		i.nchans = 1;
	end
	i.nsamples = (i.nbytes-i.bytes_to_data)/i.sample_bytes;
	i.totalsamples = i.nsamples * i.nchans;
    i.sr=1;
	return
case {'long', 'int32'};
	i.sample_bytes = 4;
	i.sample_type = 'int32';
	if isfield(i, 'bytes_to_data')
		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	if ~isfield(i, 'nchans')
		i.nchans = 1;
	end
	i.nsamples = (i.nbytes-i.bytes_to_data)/i.sample_bytes;
	i.totalsamples = i.nsamples * i.nchans;
    i.sr=1;
	return
case {'float', 'float32'};
	i.sample_bytes = 4;
	i.sample_type = 'float32';
	if isfield(i, 'bytes_to_data')
		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	if ~isfield(i, 'nchans')
		i.nchans = 1;
	end
	i.nsamples = (i.nbytes-i.bytes_to_data)/i.sample_bytes;
	i.totalsamples = i.nsamples * i.nchans;
    i.sr=1;
	return
case {'double', 'float64'};
	i.sample_bytes = 8;
	i.sample_type = 'float64';
	if isfield(i, 'bytes_to_data')
 		fseek(fd, i.bytes_to_data, 0);
	else
		i.bytes_to_data = 0;
	end
	if ~isfield(i, 'nchans')
		i.nchans = 1;
	end
	i.nsamples = (i.nbytes-i.bytes_to_data)/i.sample_bytes;
	i.totalsamples = i.nsamples * i.nchans;
    i.sr=1;
	return
otherwise
	error(['unknown format: >' i.format '<']);
end
