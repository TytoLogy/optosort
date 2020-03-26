%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% BEGIN testThreshold.m 9 Mar 2020
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% save data to mat file?
SAVEMAT = 0;


% data locations

%{
% Data for 1382_20191212_02_02_3200 has good recordings on 
% Channels 4, 5, 11, 14
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
% rawFile = '1382_20191212_02_02_3200_FRA.dat';
Channel = 4;
%}

rawPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort/IC-probe';
rawFile = '1407_20200305_01_01_550_FRA.dat';
Channel = 8;

fI = OptoFileName(fullfile(rawPath, rawFile));

%------------------------------------------------------------------------
%% get filtered data
%------------------------------------------------------------------------
% probably don't need full raw data in D....
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', Channel);
%% store working data to save time in future
if SAVEMAT
	tmpF = fullfile('~/Work/Data/TestData/MT', ...
					[fI.base sprintf('_Chan%d_TracesByStim.mat', Channel)]); %#ok<UNRCH>
	save(tmpF, ...
		'Dinf', 'tracesByStim');
end

%% test some things from threshold_opto_data
% check reps
nstim = numel(tracesByStim);
[trRows, trCols] = size(tracesByStim);

reps_by_stim = zeros(trRows, trCols);
for r = 1:trRows
	for c = 1:trCols
		[reps_by_stim(r, c), ~] = size(tracesByStim{r, c}');
% 		fprintf('tracesByStim{%d, %d} has %d reps\n', r, c, reps_by_stim(r, c));
	end
end
if length(unique(reshape(reps_by_stim, numel(reps_by_stim), 1))) > 1
	warning('%s: unequal number of reps in tracesByStim', mfilename);
	nreps = max((reshape(nreps, numel(reps_by_stim), 1)));
	fprintf('Using max nreps for allocation: %d\n', nreps);
else
	nreps = reps_by_stim(1);
end

%% calculate rms parameters for thresholding
netrmsvals = zeros(nstim, nreps);
maxvals = zeros(nstim, nreps);
% find rms, max vals for each stim
for s = 1:nstim
	netrmsvals(s, :) = rms(tracesByStim{s});
	maxvals(s, :) = max(abs(tracesByStim{s}));
end
% compute overall mean rms for threshold
fprintf('Calculating mean and max RMS for data...\n');
mean_rms = mean(reshape(netrmsvals, numel(netrmsvals), 1));
fprintf('\tMean rms: %.4f\n', mean_rms);
% find global max value (will be used for plotting)
global_max = max(max(maxvals));
fprintf('\tGlobal max abs value: %.4f\n', global_max);


cInfo = CurveInfo(Dinf);


%------------------------------------------------------------------------
%% threshold data
%------------------------------------------------------------------------
sep_print('threshold_opto_data...');
[spiketimes, cInfo, threshSettings] = threshold_opto_data(Dinf, tracesByStim);

%% compare FRA
[spiketimes, cInfo, threshSettings] = threshold_opto_data(Dinf, tracesByStim);
[spiketimesO, ~, ~] = threshold_opto_data(Dinf, tracesByStim, 'FRAMETHOD', 'OLD');

% window for spike count
frawin = [cInfo.Dinf.audio.Delay (cInfo.Dinf.audio.Delay + cInfo.Dinf.audio.Duration)];

varlist = cInfo.varlist;
FRA = computeFRA(spiketimes, varlist{1}, varlist{2}, frawin);
FRAo = computeFRA(spiketimesO, varlist{1}, varlist{2}, frawin);
hFRA = plotFRA(FRA, 'dB');
hFRAo = plotFRA(FRAo, 'dB');
 % set plot name
set(hFRA, 'Name', 'spiketimes');
set(hFRAo, 'Name', 'spiketimesO');

%% extract spike waveform "snippets"

% snippet window: [preTS time, postTStime], milliseconds
SnippetWindow = [1 2];


snips = cell(cInfo.ntrials, 1);
for tr = 1:cInfo.ntrials
	snips{tr} = extract_snippets(spiketimes{tr}, tracesByStim{tr}, ...
												SnippetWindow, cInfo.Dinf.indev.Fs);
end

%% 

figure(1);

t = 

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% END testThreshold.m
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------