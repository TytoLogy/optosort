%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
rawPath = '~/Work/Data/TestData/MT';
% rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
rawFile = '1382_20191212_02_02_3200_FRA.dat';
% Channels = [4, 5, 11, 14];

%------------------------------------------------------------------------
% get filtered data
%------------------------------------------------------------------------
% probably don't need full raw data in D....
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', 4);

%------------------------------------------------------------------------
%% threshold data
%------------------------------------------------------------------------
sep_print('threshold_opto_data...');
[spiketimes, cI, threshSettings] = threshold_opto_data(D, Dinf, tracesByStim);


%% extract spike waveform "snippets"

% snippet window: [preTS time, postTStime], milliseconds
SnippetWindow = [1 2];

snips = cell(cInfo.ntrials, 1);
for tr = 1:cInfo.ntrials
	snips{tr} = extract_snippets(spiketimes{tr}, tracesByStim{tr}, ...
												SnippetWindow, cInfo.Dinf.indev.Fs);
end