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

% S object file
Sfile = '~/Work/Data/TestData/MT/1382_20191212_02_02_3200_Sobj.mat';

% load data
load(Sfile)

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

%% get list of units
unitList = S.listUnits;

fprintf('Units in %s:\n', Sfile);
fprintf('\t%d', unitList);
fprintf('\n');

%% stimcache object

for f = 1:S.Info.nFiles
	Dinf = S.Info.FileData(f).Dinf;
	
	CI{f} = CurveInfo(Dinf);
% 	switch upper(Dinf.test.Type)
% 		case 'FREQ'
% 			% frequency tuning curve
% 				SC{f} = FREQcache(S.Info.getStimulusCacheByFile(f), ...
% 										Dinf.test.Type, Dinf.test.Name);
% 		case 'LEVEL'
% 				SC{f} = LEVELcache(S.Info.getStimulusCacheByFile(f), ...
% 										Dinf.test.Type, Dinf.test.Name);
% 		otherwise
% 			error('unsupported test type %s', Dinf.test.Type);
% 	end

end