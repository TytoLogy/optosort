%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
% rawFile = '1382_20191212_02_02_3200_FRA.dat';
% Channels = [4, 5, 11, 14];

rawPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort/IC-probe';
rawFile = '1407_20200305_01_01_550_FRA.dat';
Channels = 8;

fI = OptoFileName(fullfile(rawPath, rawFile));

%------------------------------------------------------------------------
%% get filtered data
%------------------------------------------------------------------------
ch = Channels(1);

% probably don't need full raw data in D....
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', ch);
%% store working data to save time in future
if SAVEMAT
	tmpF = fullfile('~/Work/Data/TestData/MT', ...
					[fI.base sprintf('_Chan%d_TracesByStim.mat', ch)]);
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
[varlist, nvars] = get_opto_varlist(cInfo)

%------------------------------------------------------------------------
%% threshold data
%------------------------------------------------------------------------
sep_print('threshold_opto_data...');
[spiketimes, cInfo, threshSettings] = threshold_opto_data(Dinf, tracesByStim);


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