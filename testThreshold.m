%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------

% save data to mat file?
SAVEMAT = 0;


% data locations

%
% Data for 1382_20191212_02_02_3200 has good recordings on 
% Channels 4, 5, 11, 14
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
Channel = 4;
%

%{
% data for FRA - this is correct data to use for testing as it was
collected using updated FRA routine in opto program

rawPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort/IC-probe';
rawFile = '1407_20200305_01_01_550_FRA.dat';
Channel = 8;
%}

fI = OptoFileName(fullfile(rawPath, rawFile));

%------------------------------------------------------------------------
%% get filtered data
%------------------------------------------------------------------------
% probably don't need full raw data in D....
nChannels = length(Channel);
Traces = cell(nChannels, 1);
for c = 1:nChannels
% 	[D, Dinf, Traces{c}] = getFilteredOptoData( ...
% 											fullfile(rawPath, rawFile), ...
% 											'Filter', [300 5000], ...
% 											'Channel', Channel(c));
	[Traces{c}, Dinf] = getOptoTracesByStim( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', Channel(c));
end
%% store working data to save time in future
if SAVEMAT
	tmpF = fullfile('~/Work/Data/TestData/MT', ...
					[fI.base sprintf('_Chan%d_TracesByStim.mat', Channel)]); %#ok<UNRCH>
	save(tmpF, ...
		'Dinf', 'tracesByStim');
end

%------------------------------------------------------------------------
%% threshold data
%------------------------------------------------------------------------
sep_print('threshold_opto_data...');
[spikedata, cInfo, threshSettings] = threshold_opto_data(Dinf, tracesByStim);
