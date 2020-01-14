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
for f = 1:nFiles
	fprintf('DataFile{%d} = %s\n', f, DataFile{f});
end
fprintf('Animal: %s\n', F(1).animal);

% !!!!!!!
% in future, will need to have user specify all the files for this
% recording session and recording location!!!
% ¡¡¡¡¡¡¡

% channel(s) of data to obtain
% channel = 8;
Channels = [11 9 14];
nChannels = length(Channels);

% filter info
BPfilt.Fc = [300 4000];
BPfilt.forder = 5;
BPfilt.ramp = 1;

%% setup for channel data conversion



%% loop through files

% allocate some things
fileIDstartbin = zeros(nFiles, 1);
fileIDendbin = zeros(nFiles, 1);
tmpFs = zeros(nFiles, 1);

for f = 1:nFiles

	% get data for each file and channel and convert to row vector format
	% algorithm:
	%		(1) put each sweep for this channel in a {1, # sweeps} cell array
	%				cSweeps
	%		(2) make a note of the length of each sweep to use for
	%				markers/timestamps
	%		(3) after cSweeps is built, convert to a row vector using cell2mat
	% 

	% read in data
	% use readOptoData to read in raw data. 
	[D, Dinf] = readOptoData(fullfile(DataPath, DataFile{f}));
	% Fix test info
	Dinf = correctTestType(Dinf);
	% build filter for neural data
	BPfilt.Fs = Dinf.indev.Fs;
	BPfilt.Fnyq = Dinf.indev.Fs / 2;
	BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
	[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');

	% check to make sure consistent # of sweeps (aka trials)
	if Dinf.nread ~= length(D)
		error('%s: mismatch in Dinf.nread (%d) and length(D) (%d)', ...
					mfilename, Dinf.nread, length(D));
	end
	
	% build into sweeps by channel format
	fprintf('Test type: %s\n', Dinf.test.Type);
	[fData(f).cSweeps, fData(f).startSweepBin, fData(f).endSweepBin] = ...
					buildChannelData(Channels, BPfilt, D, Dinf); %#ok<*SAGROW>
	% check the start and end sweep bin data for consistency
	if check_sweeps(fData(f).startSweepBin)
		warning('Inconsistent startSweepBin values across channels!!!!');
	end
	if check_sweeps(fData(f).endSweepBin)
		warning('Inconsistent endSweepBin values across channels!!!!');
	end
	% store sample for start of this file (should be 1); use channel 1 value
	fData(f).fileStartBin = fData(f).startSweepBin{1}(1);
	% store sample for end of this file
	fData(f).fileEndBin = fData(f).endSweepBin{1}(end);
	% store file information struct
	fData(f).Dinf = Dinf;
	% to avoid any issues, should make sure sample rates are consistent
	% need in implement a check somehow.  this code doesn't work:
	% if any(Dinf.indev.Fs ~= [fData(:).Dinf.indev.Fs])
	%	error('Sample Rate Mismatch found!');
	% end
	tmpFs(f) = fData(f).Dinf.indev.Fs;
	
	% calculate start and end bins for each file's data
	if f == 1
		fileStartBin(f) = fData(f).fileStartBin;
		fileEndBin(f) = fData(f).fileEndBin;
	else
		% add 1 to prior end bin for start
		fileStartBin(f) = fileStartBin(f-1) + 1;
		fileEndBin(f) = fileStartBin(f) + fData(f).fileEndBin - 1;
	end
end

% check sampling rates
if ~all(tmpFs(1) == tmpFs)
	error('Sample Rate Mismatch!!!');
else
	% store overall sample rate
	Fs = tmpFs(1);
end

%% convert cSweeps to vector (or matrix) for each channel
% to save on memory requirements, clear channel data after adding to nex
% data struct.

% create output .nex file name
NexFileName = [F.base '.nex'];
% start new nex file data
nD = nexCreateFileData(Fs);


% channel
c = 1;
% this will be a matrix of format
% 	[# channels, (# sweeps) * (# samples per sweep)
cVector = cell(1, nFiles);
for f = 1:nFiles
	cVector{1, f} = fData(f).cSweeps(c, :);
end
% concatenate cell array and convert to vector
% steps:
%	concatenate: tmp = [cVector{:}];
%	tmpVector = cell2mat(tmp);

nD = nexAddContinuous(nD, startSweepTime{c}(1), Dinf.indev.Fs, ...
									cell2mat([cVector{:}]), sprintf('spikechan_%d', Channels(c)));


%% Convert sweep samples to time (seconds)


for c = 1:nChannels
	cVector{c} = cell2mat(cSweeps(c, :));
	% convert samples to timestamps; need to subtract 1 to start at time = 0
	startSweepTime{c} = (startSweepBin{c} - 1) / Dinf.indev.Fs;
	endSweepTime{c} = (endSweepBin{c} - 1) / Dinf.indev.Fs;
end

%% do a check on sweep bins
if check_sweeps(startSweepBin)
	warning('Inconsistent startSweepBin values across channels!!!!');
else
	fprintf('startSweepBin values are consistent across channels\n');
end
if check_sweeps(endSweepBin)
	warning('Inconsistent endSweepBin values across channels!!!!');
else
	fprintf('endSweepBin values are consistent across channels\n');
end
%% write output
% create output .nex file name
NexFileName = [F.base '.nex'];
% start new nex file data
nD = nexCreateFileData(Dinf.indev.Fs);

% add each channel's data as a continuous variable
%{
specify start time (t(1)), digitizing frequency (Fs), data and name
[nexFile] = nexAddContinuous( nexFile, startTime, adFreq, values, name ) 
        -- adds continuous variable to nexFile data structure

INPUT:
  nexFile - nex file data structure created in nexCreateFileData
  startTime - time of the first data point in seconds
  adFreq - A/D sampling rate of continuous variable in samples per second
  values - vector of continuous variable values in milliVolts
  name - continuous variable name  
%}
for c = 1:nChannels
	nD = nexAddContinuous(nD, startSweepTime{c}(1), Dinf.indev.Fs, ...
									cVector{c}, sprintf('spikechan_%d', Channels(c)));
end								
% add start sweep time stamps as event - assume consistent across channels!
%{
[nexFile] = nexAddEvent( nexFile, timestamps, name ) -- adds an event 
            (series of timestamps) to nexFile data structure

INPUT:
  nexFile - nex file data structure created in nexCreateFileData
  timestamps - vector of event timestamps in seconds
  name - event name
%}
% events must be in column format...?
nD = nexAddEvent(nD, startSweepTime{1}, 'startsweep');
% add end sweep time stamps as event - assume consistent across channels!
nD = nexAddEvent(nD, endSweepTime{1}, 'endsweep');

% write to nexfile
writeNexFile(nD, NexFileName);


