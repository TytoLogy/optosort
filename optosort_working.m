%------------------------------------------------------------------------
% optosort_working.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% working script for opto data <--> plexon
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 8 January, 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%% Initial things

sepstr = '----------------------------------------------------';
fprintf('\n%s\n', sepstr);
fprintf('%s running...\n', mfilename);
fprintf('%s\n', sepstr);

% add path to .nex file utils
if ~exist('nexCreateFileData', 'file')
	fprintf('Adding NEX file utility paths\n')
	addpath(['~/Work/Code/Matlab/stable/Toolbox/NeuroExplorer' ...
					'/HowToReadAndWriteNexAndNex5FilesInMatlab']);
end
% add path to opto if needed
if ~exist('readOptoData', 'file')
	fprintf('Adding opto paths\n')
	addOptoPaths
end

% define path to data file and data file for testing
defineSampleData
fprintf('DataPath = %s\n', DataPath);
fprintf('DataFile = %s\n', DataFile);
fprintf('Animal: %s\n', F.animal);

% !!!!!!!
% in future, will need to have user specify all the files for this
% recording session and recording location!!!
% ¡¡¡¡¡¡¡

% channel(s) of data to obtain
% channel = 8;
Channels = [11 9 14];

% filter info
BPfilt.Fc = [100 5000];
BPfilt.forder = 3;
BPfilt.ramp = 1;

%% read in data

% can probably just use readOptoData
% % the data will be further processed during sorting. for now, use a fairly
% % broad filter
% filtband = [5000 10000];
% % get the data - only need the raw data (stored in D cell array of each
% % stimulus presentation) and the information in Dinf about the test.
% [D, Dinf] = getFilteredOptoData(fullfile(datadir, datafile), ...
% 											'filter', filtband, ...
% 											'channel', channel);
% use readOptoData to read in raw data. 
[D, Dinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
Dinf = correctTestType(Dinf);
fprintf('Test type: %s\n', Dinf.test.Type);

%% quick plot of data for one trial

% need colorspace to distinguish channels
lineStyles = linspecer(16);
figure(10)
axes('NextPlot','replacechildren', 'ColorOrder', lineStyles);
plot(D{1}.datatrace);
tl = cell(1, 16);
for c = 1:16
	tl{c} = num2str(c);
end
legend(tl);

%% setup output

% build filter for neural data
BPfilt.Fs = Dinf.indev.Fs;
BPfilt.Fnyq = Dinf.indev.Fs / 2;
BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');


% create output .nex file name
NexFile = [F.base '.nex'];
% start new nex file data
nD = nexCreateFileData(Dinf.indev.Fs);

% check to make sure consistent # of sweeps (aka trials)
if Dinf.nread ~= length(D)
	error('%s: mismatch in Dinf.nread (%d) and length(D) (%d)', ...
				mfilename, Dinf.nread, length(D));
end

% get data for each channel
% algorithm:
%		(1) put each sweep for this channel in a {1, # sweeps} cell array
%				cSweeps
%		(2) make a note of the length of each sweep to use for
%				markers/timestamps
%		(3) after cSweeps is built, convert to a row vector using cell2mat
% initialize sweep indices to empty arrays

[cSweeps2, s1, s2] = buildChannelData(Channels, BPfilt, D, Dinf);


%%
sweepStartSample = [];
sweepEndSample = [];
% loop through channels
for c = 1:length(Channels)
	channel = Channels(c);
	% initialize cSweeps to a row cell "vector"
	cSweeps = cell(1, Dinf.nread);
	% loop through each sweep
	for s = 1:Dinf.nread
		% assign channel sweep data to cSweeps after filtering
		% need to transpose to row vector
		cSweeps{1, s} = D{s}.datatrace(:, channel)';
		% remove initial offset
		cSweeps{1, s} = cSweeps{1, s} - cSweeps{1, s}(1);
		% filter data
		cSweeps{1, s} = filtfilt(BPfilt.b, BPfilt.a, ...
									sin2array(cSweeps{1, s}, ...
										Dinf.indev.Fs, BPfilt.ramp));
		% plot raw and filtered data
		plot(D{s}.datatrace(:, channel)', 'k');
		hold on
			plot(cSweeps{1, s}, 'b');
		hold off
		drawnow
	end

	% build list of sweep start and end indices (in units of samples)
	if isempty(sweepStartSample)
		% initialize sweepStartSample and sweepEndSample to store stop/start
		% locations
		sweepStartSample = zeros(1, Dinf.nread);
		sweepEndSample = zeros(1, Dinf.nread);
		% store index points
		if s ~= 1
			sweepStartSample(s) = sweepEndSample(s-1) + 1;
		else
			sweepStartSample(s) = 1;
		end
		sweepEndSample(s) = sweepStartSample(s) + length(cSweeps{1, s});
	end
end


