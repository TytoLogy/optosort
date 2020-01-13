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
nChannels = length(Channels);

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

%% setup for channel data conversion

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

%% get data for each channel and convert to row vector format
% algorithm:
%		(1) put each sweep for this channel in a {1, # sweeps} cell array
%				cSweeps
%		(2) make a note of the length of each sweep to use for
%				markers/timestamps
%		(3) after cSweeps is built, convert to a row vector using cell2mat
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This will have to be modified to deal with multiple files!!!!
% ¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡
[cSweeps, startSweepBin, endSweepBin] = buildChannelData(Channels, BPfilt, ...
																				D, Dinf);
% convert cSweeps to vector (or matrix)
% this will be a matrix of format
% 	[# channels, (# sweeps) * (# samples per sweep)
cVector = cell(nChannels, 1);
startSweepTime = cell(nChannels, 1);
endSweepTime = cell(nChannels, 1);
for c = 1:nChannels
	cVector{c} = cell2mat(cSweeps(c, :));
	% convert samples to timestamps; need to subtract 1 to start at time = 0
	startSweepTime{c} = (startSweepBin{c} - 1) / Dinf.indev.Fs;
	endSweepTime{c} = (endSweepBin{c} - 1) / Dinf.indev.Fs;
end

%% do a check on start sweep times
tmpStart = cell2mat(startSweepBin);
tmpEnd = cell2mat(endSweepBin);

if sum((sum(tmpStart - tmpStart(1, :))))
	warning('Inconsistent startSweepBin values across channels!!!!');
else
	fprintf('startSweepBin values are consistent across channels\n');
end
if sum((sum(tmpEnd - tmpEnd(1, :))))
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
nD = nexAddContinuous(nD, startSweepTime{c}(1), Dinf.indev.Fs, ...
									cVector{c}, sprintf('spikechan_%d', Channels(c)));
								
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


