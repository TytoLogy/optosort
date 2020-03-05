%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
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


%% test snippetizing


%% load test data
load snipindata.mat
Fs = cInfo.Dinf.indev.Fs;
timestamps = spiketimes{1};
sweeps = tracesByStim{1};

% snippet window: [preTS time, postTStime], milliseconds
SnippetWindow = [1 2];


snips = extract_snippets(timestamps, sweeps, SnippetWindow, Fs)