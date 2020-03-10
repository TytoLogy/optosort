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
Channels = [4 5];
%

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

%%
tmpF = fullfile('~/Work/Data/TestData/MT', ...
				fI.newname(  sprintf('Chan%d-%d_TracesByStim', ...
															Channels(1), Channels(end)), ...
								'mat') );
%% store working data to save time in future
if SAVEMAT
	tmpF = fullfile('~/Work/Data/TestData/MT', ...
				fI.newname(  sprintf('Chan%d-%d_TracesByStim', ...
															Channels(1), Channels(end)), ...
								'mat') ); %#ok<UNRCH>
	save(tmpF, 'Dinf', 'Traces', 'cInfo');
end

%------------------------------------------------------------------------
%% threshold data
%------------------------------------------------------------------------
sep_print('threshold_opto_data...');
spikedata = cell(nChannels, 1);
tset = cell(nChannels, 1);
for c = 1:nChannels
	fprintf('Thresholding channel %d\n', Channels(c));
	[spikedata{c}, tset{c}] = threshold_opto_data(cInfo, Traces{c});
end

%% show thresholded data

spikedata{1}.ts{1}
spikedata{1}.snips{1}



