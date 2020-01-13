function [cD, startI, endI] = buildChannelData(Channels, BPfilt, D, Dinf, varargin)
%------------------------------------------------------------------------
% Output = function_template(Input)
%------------------------------------------------------------------------
% <project>:<subproject>:<function_name>
%------------------------------------------------------------------------
% 
% get data for each channel
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Channels		list of A/D (spike) channels to export
%	BPfilt		struct specifying bandpass filter settings
% 					if BPfilt is empty ([]), no filtering will be performed
% 			required fields:
% 			forder	filter order
% 			ramp		onset/offset signal ramp in milliseconds
% 			Fs			signal sampling rate in samples/sec
% 			b			bandpass filter coefficients (e.g., from butter())
% 			a			bandpass filter coefficients (e.g., from butter())
% 	D				raw data from opto program output file
% 	Dinf			data information struct from opto program output file
% 	
% 	Optional Input Args:
% 		'REMOVE_OFFSET'	<'Y'/'N'>	remove DC offset from sweeps?
% 		'PLOT_SWEEPS'		<'N'/'Y'>	plot sweeps?
%
% Output Arguments:
% 	cD			{# Channels, # sweeps} cell array of channel sweep data
% 				each element of cD will be a row vector of raw data for a sweep
% 	startI	[# sweeps] vector holding start sample timestamp for each sweep
% 	endI		[# sweeps] vector holding end samples timestamp for each sweep
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 January, 2020 (SJS)
%
% Revisions:
%	13 Jan 2020 (SJS): added comments , optional args
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% settings
%------------------------------------------------------------------------
deOffset = true;
plotSweeps = false;
filterData = true;

if ~isempty(varargin)
	argIndx = 1;
	while argIndx <= length(varargin)
		fprintf('%s\n', upper(varargin{argIndx}))
		switch upper(varargin{argIndx})
			case {'REMOVE_OFFSET'}
				tmpArg = varargin{argIndx + 1};
				if strcmpi(tmpArg(1), 'Y')
					deOffset = true;
				else
					deOffset = false;
				end
				argIndx = argIndx + 2;
			case {'PLOT_SWEEPS'}
				tmpArg = varargin{argIndx + 1};
				if strcmpi(tmpArg(1), 'Y')
					plotSweeps = true;
				else
					plotSweeps = false;
				end
				argIndx = argIndx + 2;
			otherwise
				error('%s: unknown input arg %s', mfilename, varargin{argIndx});
		end
	end
end
	
if isempty(BPfilt)
	filterData = false;
end

%------------------------------------------------------------------------
% process data
%------------------------------------------------------------------------

% initialize sweep indices to empty arrays
startI = [];
endI = [];

% initialize cD to a store sweeps for each channel
cD = cell(length(Channels), Dinf.nread);

% loop through channels
for c = 1:length(Channels)
	channel = Channels(c);
	% loop through each sweep
	for s = 1:Dinf.nread
		% assign channel sweep data to cD after filtering
		% need to transpose to row vector
		cD{c, s} = D{s}.datatrace(:, channel)';
		if deOffset
			% remove initial offset
			cD{c, s} = cD{c, s} - cD{c, s}(1);
		end
		
		% %%%% DEBUG
		cD{c, s}(1) = exp(1);
		cD{c, s}(end) = -exp(1);
		
		% filter data
		if filterData
			cD{c, s} = filtfilt(BPfilt.b, BPfilt.a, ...
										sin2array(cD{c, s}, ...
										Dinf.indev.Fs, BPfilt.ramp));
		end
		% plot raw and filtered data
		if plotSweeps
			plot(D{s}.datatrace(:, channel)', 'k');
			hold on
				plot(cD{c, s}, 'b');
			hold off
			title(sprintf('Channel: %d  Sweep: %d(%d)', channel, s, Dinf.nread));
			drawnow
		end
	end

	% build list of sweep start and end indices (in units of samples)
	if isempty(startI)
		% initialize startI and endI to store stop/start
		% locations
		startI = zeros(1, Dinf.nread);
		endI = zeros(1, Dinf.nread);
		% store index points
		if s ~= 1
			startI(s) = endI(s-1) + 1;
		else
			startI(s) = 1;
		end
		endI(s) = startI(s) + length(cD{c, s});
	end
end


