% get spikes by stimulus
% script to work out sorting spikes by stimulus

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 19 February, 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%--------------------------------------------------------------------------

%{
how to do this...

one way:
	separate spikes by file (in other words, by test type - freq resp, rate
	level, wav stim etc)
	then, within each group, match up stimulus type from stimcache with the
	appropriate spikes for each "sweep"

other way:
	no need to separate by file, 
	merge the stimcache arrays to generate a master merged list of stimulus
	by sweep
	then, sort the spikes into spikes by sweep (for each channel)
	can then search through stimcache to get appropriate indices for desired
	stimuli. 
	or separate spikes by stimulus... but this will disrupt sweep-based
	analysis....

	idea: could create a table form of stimcache to ease searching/retrieval
	of stimulus types?

problem: each test has different parameters: e.g., freq-tuning varies
freq, freq-level varies level, FRA both freq and level, WAV stimulus and
level.... how to have a common basis for this to put into uniform table?
also... OPTO!!!????!!

maybe best to separate by file, keep spiketimes referenced as in MERGED
file and then have separate processing stages depending on test type...
%}

%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
sortedPath = '~/Work/Data/TestData/MT';
rawPath = '~/Work/Data/TestData/MT';
nexPath = '~/Work/Data/TestData/MT';

% sorted data file
% sortedFile = '1323_20190722_03_02_632_MERGE.mat';
sortedFile = '1382_20191212_02_02_3200.mat';

% nexinfo file
% nexInfoFile = '1323_20190722_03_02_632_MERGE_nexinfo.mat';
nexInfoFile = '1382_20191212_02_02_3200_MERGE_nexinfo.mat';

% nex file
nexFile = '1382_20191212_02_02_3200_MERGE.nex';

S = load_plexon_data(fullfile(sortedPath, sortedFile), ...
							fullfile(sortedPath, nexInfoFile));

%% get stimcache data

% these are located in the S.Info.FileData(1).Dinf.test.stimcache
stimcache = cell(S.Info.nFiles, 1);
for f = 1:S.Info.nFiles
	stimcache{f} = S.Info.getStimulusCacheByFile(f);
end

%% Get Spikes for this file
f = 1;

% get spikes by sweep. this will return a cell array of tables with size
% {nSweeps, 1}. Spike times in the table will be aligned to sweep onset
fS = S.spikesForAnalysis(f, 'sweep');
% get spikes in orginal form (aligned to start of merged file)
fS_orig = S.spikesForAnalysis(f);

%% stimcache object

for f = 1:S.Info.nFiles
	Dinf = S.Info.FileData(f).Dinf;
	
	CI{f} = CurveInfo(Dinf);
end

%% sweep data by stimulus

f = 1;
[stimIndex, stimVar] = CI{f}.getStimulusIndices;
fprintf('Found spike data for units: ')
fprintf('%d  ', S.listUnits);
fprintf('\n\n');

allUnits = S.listUnits;
validUnits = allUnits(allUnits > 0);
spikesByStim = fS{stimIndex{1}};
curveInfo = CI{f};

% curvedata by test type
u = 1;
unit = validUnits(u);
switch CI{f}.testtype
	case 'FREQ'
		
		
	otherwise
		error('%s: curve %s not implemented', mfilename, CI{f}.testtype)
end

