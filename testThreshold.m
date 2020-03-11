%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------

% save data to mat file?
SAVEMAT = 0;

% spike window (pre, post timestamp, in milliseconds)
SpikeWindow = [1 2];

% data locations

%
% Data for 1382_20191212_02_02_3200 has good recordings on 
% Channels 4, 5, 11, 14
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
Channels = [4 5];
%}

%{
% Test WAV data
% 1382_20191212_02_02_3200_WAV.dat
% 1382_20191212_02_02_3200_WAV_PSTH.fig
% 1382_20191212_02_02_3200_WAV_wavinfo.mat
% Channels 4, 5, 11, 14
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_WAV.dat';
Channels = [4 5];
%}

%{
% data for FRA - this is correct data to use for testing as it was
collected using updated FRA routine in opto program

rawPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort/IC-probe';
rawFile = '1407_20200305_01_01_550_FRA.dat';
Channels = 8;
%}

% create opto file name object (helps to create new file names for these
% data)
fI = OptoFileName(fullfile(rawPath, rawFile));

%------------------------------------------------------------------------
%% get filtered data
%------------------------------------------------------------------------
% probably don't need full raw data in D....
nChannels = length(Channels);
Traces = cell(nChannels, 1);
for c = 1:nChannels
% 	[D, Dinf, Traces{c}] = getFilteredOptoData( ...
% 											fullfile(rawPath, rawFile), ...
% 											'Filter', [300 5000], ...
% 											'Channels, Channels(c));
	sep_print(sprintf('testThreshold: Reading data for A/D channel %d', ...
																			Channels(c)));
	if c == 1	
		[Traces{c}, Dinf] = getOptoTracesByStim( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', Channels(c));
		% create CurveInfo struct
		cInfo = CurveInfo(Dinf);
	else
		Traces{c} = getOptoTracesByStim( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', Channels(c));
	end
end

%------------------------------------------------------------------------
%% threshold data and plot detected waveforms
%------------------------------------------------------------------------
sep_print('threshold_opto_data...');
spikes = cell(nChannels, 1);
tset = cell(nChannels, 1);
for c = 1:nChannels
	fprintf('Thresholding channel %d\n', Channels(c));
	[spikes{c}] = threshold_opto_data(cInfo, Traces{c}, ...
													'Method', 'RMS', ...
													'Threshold', 5, ...
													'Spike_Window', SpikeWindow);
end

%  plot waveforms
for c = 1:nChannels
	figure(c)
	SD = spikes{c};

	plot_snippets(SD, SpikeWindow, cInfo.ADFs);
	title({cInfo.F.file, sprintf('Channel %d', Channels(c))}, ...
							'Interpreter', 'none');


	set(gcf, 'Name', [fI.base '-Ch' num2str(Channels(c))])
end

%% store working data to save time in future
if SAVEMAT
	tmpF = fullfile('~/Work/Data/TestData/MT', ...
				fI.newname(  sprintf('Chan%d-%d_TracesByStim', ...
															Channels(1), Channels(end)), ...
								'mat') ); %#ok<UNRCH>
	save(tmpF, 'Dinf', 'Traces', 'cInfo', 'Channels', 'spikes', 'tset');
end

%% plot data in explorable form

tracesByStim = Traces{1};

f = figure('Units','Normalized','Position',[0.25 0.25 0.5 0.5]);
a =   axes('Units','Normalized','Position',[0.05 0.15, 0.75 0.75]);




