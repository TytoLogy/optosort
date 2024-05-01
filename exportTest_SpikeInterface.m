%------------------------------------------------------------------------
% exportTest_SpikeInterface.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% sample script to show how to export data for spike sorting using 
% SpikeInterface (and Kilosort?)
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 1 May, 2024 from exportTestRaw (SJS)
%
% Revisions:
%  14 May, 2020 exportTestRaw created from exportTest (SJS)
%   1 May, 2024 exportTest_SpikeInterface created from exportTestRaw (SJS)
%               uses OutputShape->'SamplesChannels' option to write output
%               data in format [samples, channels] needed for Kilosort
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% check to make sure tytoLogy (esp opto) things are on path
if ~exist('readOptoData', 'file')
	fprintf(['readOptoData (and possibly other TytoLogy/opto functions' ...
						'not found!\n'])
	fprintf('Please check Matlab path(s)\n');
	fprintf(['e.g.:\n' ...
				'addpath(''~/Work/Code/Matlab/dev/TytoLogy/Experiments/Opto'')\n']);
end

if ~exist('SpikeData', 'file')
	fprintf('SpikeData.m class definition file not found!\n')
	fprintf('This is usually found in the OptoObjects folder\n')
	fprintf('Please check Matlab path(s)\n');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% PATHS TO DATA FILES
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% can specify individual file paths... these must match the paths in
% exportOpts.DataFile!!!!!!
% exportOpts.DataPath = {'~/Work/Data/TestData/MT_IC';
% 				'~/Work/Data/TestData/MT_IC'; ...
% 				'~/Work/Data/TestData/MT_IC'; };
% or single path:
% exportOpts.DataPath = '~/Work/Data/TestData/MT_IC';
% exportOpts.DataPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort';
% exportOpts.DataPath = '/Volumes/SJS_XFER/Work/MT-IC-R-data';
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% DATA FILES
%------------------------------------------------------------------------

%---------------------------------------
% IC data from probe, include WAV test
% exportOpts.DataPath = ['/Volumes/Wenstrup Laboratory/By User/SJS/Data' ...
% 							'/SpikeSort/IC-probe/1372'];
% exportOpts.DataFile = {'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
% 				'1372_20191126_03_01_1500_BBN.dat'; ...
% 				'1372_20191126_03_01_1500_FRA.dat'; ...
% 				'1372_20191126_03_01_1500_WAV.dat'; };
% exportOpts.TestFile = {'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
% 				'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
% 				'1372_20191126_03_01_1500_FRA_testdata.mat'; ...
% 				''; };
% exportOpts.Channels = [11 9 14];
%---------------------------------------
% IC data from probe, NO WAV test
% exportOpts.DataPath = ['/Volumes/Wenstrup Laboratory/By User/SJS/Data' ...
% 							'/SpikeSort/IC-probe/1372'];
% exportOpts.DataPath = '~/Work/Data/TestData/MT';
% exportOpts.DataFile = {	'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
% 								'1372_20191126_03_01_1500_BBN.dat'; ...
% 								'1372_20191126_03_01_1500_FRA.dat'; };
% exportOpts.TestFile = {	'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
% 								'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
% 								'1372_20191126_03_01_1500_FRA_testdata.mat'; };
% exportOpts.Channels = [11 9 14];
% you can specify an output path and nex file name, or just leave blank
% and export_plexon_data will create one in current directory
% exportOpts.OutputPath = exportOpts.DataPath;
% exportOpts.OutputFile = '1372_20191126_03_01_1500_test.nex';
%---------------------------------------

%---------------------------------------
% MG sample data with probe
% 1382_20191212_02_02_3200_WAV.dat
% 1382_20191212_02_02_3200_WAV_PSTH.fig
% 1382_20191212_02_02_3200_WAV_wavinfo.mat
% exportOpts.DataPath = '~/Work/Data/TestData/MT';
% exportOpts.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'; ...
% 					'1382_20191212_02_02_3200_BBN.dat'; ...
% 					'1382_20191212_02_02_3200_FRA.dat';	};
% exportOpts.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'; ...
% 					'1382_20191212_02_02_3200_BBN_testdata.mat'; ...
% 					'1382_20191212_02_02_3200_FRA_testdata.mat'; };
% exportOpts.Channels = [4, 5, 11, 14];
% you can specify an output path and nex file name, or just leave blank
% and export_plexon_data will create one in current directory
% exportOpts.OutputPath = exportOpts.DataPath;
% exportOpts.OutputFile = '1382_20191212_02_02_3200_test.nex';
%---------------------------------------

%---------------------------------------
% IC sample data tungsten (single channel)
% 1323_20190722_03_02_632_WAV.dat
% 1323_20190722_03_02_632_WAV_testdata.mat
% 1323_20190722_03_02_632_WAV_wavinfo.mat
% exportOpts.DataPath = ['/Volumes/Wenstrup Laboratory/By User/SJS/Data' ...
% 						'/SpikeSort/IC-tungsten/'];
% exportOpts.DataFile = {	'1323_20190722_03_02_632_FREQ_TUNING.dat'; ...
% 					'1323_20190722_03_02_632_BBN.dat'; ...
% 					'1323_20190722_03_02_632_TONE_LEVEL.dat';	};
% exportOpts.TestFile = {	'1323_20190722_03_02_632_FREQ_TUNING_testdata.mat'; ...
% 					'1323_20190722_03_02_632_BBN_testdata.mat'; ...
% 					'1323_20190722_03_02_632_TONE_LEVEL_testdata.mat'; };
% exportOpts.Channels = [8];
%---------------------------------------

%---------------------------------------
% testing object
%---------------------------------------
% MG sample data with probe
% 1382_20191212_02_02_3200_WAV.dat
% 1382_20191212_02_02_3200_WAV_PSTH.fig
% 1382_20191212_02_02_3200_WAV_wavinfo.mat
% exportOpts.DataPath = '~/Work/Data/TestData/MT';
% exportOpts.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'; ...
% 					'1382_20191212_02_02_3200_BBN.dat'; ...
% 					'1382_20191212_02_02_3200_FRA.dat';	};
% exportOpts.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'; ...
% 					'1382_20191212_02_02_3200_BBN_testdata.mat'; ...
% 					'1382_20191212_02_02_3200_FRA_testdata.mat'; };
% exportOpts.Channels = [4, 5, 11, 14];
% exportOpts.OutputPath = exportOpts.DataPath;
% exportOpts.OutputFile = '1382_20191212_02_02_3200_MERGE.nex';
%---------------------------------------

%---------------------------------------
%---------------------------------------
% Example to export a single test's data
%---------------------------------------
% exportOpts.DataPath = '~/Work/Data/TestData/MT';
% exportOpts.OutputPath = exportOpts.DataPath;
% % exportOpts.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'	};
% % exportOpts.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'};
% % exportOpts.OutputFile = '1382_20191212_02_02_3200_FREQ_TUNING.nex';
% exportOpts.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'	};
% exportOpts.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'};
% exportOpts.OutputFile = '1382_20191212_02_02_3200_FREQ_TUNING.nex';
% exportOpts.Channels = [4, 5, 11, 14];
%---------------------------------------
%---------------------------------------


%---------------------------------------
%---------------------------------------
% testing scenarios
%---------------------------------------
%{
exportOpts.DataPath = '~/Work/Data/TestData/MT';
exportOpts.OutputPath = '/Volumes/ACEData/TestData/NEXtest';
exportOpts.Channels = [4, 5, 11, 14];
%}

%{
% WAV, FREQ merged
exportOpts.DataFile = {	'1382_20191212_02_02_3200_WAV.dat', ...
								'1382_20191212_02_02_3200_FREQ_TUNING.dat'};
exportOpts.TestFile = {	'1382_20191212_02_02_3200_WAV_testdata.dat', ...
								'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'};
exportOpts.OutputFile = 'FREQWAVMERGE.nex';
%}

%{
% FREQ only
exportOpts.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'};
exportOpts.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'};
exportOpts.OutputFile = 'FREQ.nex';
%}

%{
% Wav Only
exportOpts.DataFile = {	'1382_20191212_02_02_3200_WAV.dat'};
exportOpts.TestFile = {	'1382_20191212_02_02_3200_WAV_testdata.mat'};
exportOpts.OutputFile = 'WAV.nex';
%}

%{
% BBN only
exportOpts.DataFile = {	'1382_20191212_02_02_3200_BBN.dat'};
exportOpts.TestFile = {	'1382_20191212_02_02_3200_BBN_testdata.mat'};
exportOpts.OutputFile = 'BBN.nex';
%}

%{
% FRA only
exportOpts.DataFile = {	'1382_20191212_02_02_3200_FRA.dat'};
exportOpts.TestFile = {	'1382_20191212_02_02_3200_FRA_testdata.mat'};
exportOpts.OutputFile = 'FRA.nex';
%}

%{
% 1407 multichannel
exportOpts.DataPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
exportOpts.DataFile = '1407_20200309_03_01_1350_BBN.dat';
exportOpts.OutputFile = '1407_20200309_03_01_1350_BBN.bin';
exportOpts.TestFile = '1407_20200309_03_01_1350_BBN_testdata.mat';
exportOpts.Channels = [4, 5, 7, 15];
%}



% Wav Only
exportOpts.DataPath = '/media/Data/NeuroData/Mouse/Raw/BLA/1500/20230213';
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/exports/SI';
exportOpts.DataFile = {	'1500_20230213_01_0_3352_WAV.dat'};
exportOpts.TestFile = {	'1500_20230213_01_0_3352_WAV_testdata.mat'};
exportOpts.OutputFile = '1500_20230213_01_0_3352_WAV.bin';
exportOpts.Channels = 1:16;
%---------------------------------------
%---------------------------------------


%------------------------------------------------------------------------
% filter parameters for raw neural data
%------------------------------------------------------------------------
%{
% [highpass lowpass] cutoff frequencies in Hz
exportOpts.BPfilt.Fc = [300 4000];
% order of filter. note that the filtfilt() function in MATLAB is used,
% so the effective order is doubled. typically use 5
exportOpts.BPfilt.forder = 5;
% ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% transients from filtering
exportOpts.BPfilt.ramp = 1;
% filter type. 'bessel' or 'butter' (for butterworth). typically use butter
% as bessel is low-pass only and we need to remove low frequency crap from
% the raw data
exportOpts.BPfilt.type = 'butter';
%}
% no filtering
exportOpts.BPfilt = [];

%------------------------------------------------------------------------
% Specify output shape as [samples, channels] as needed for SpikeInterface
% i.e., samples in rows, channels in columns
%------------------------------------------------------------------------
exportOpts.OutputShape = 'SamplesChannels';

%------------------------------------------------------------------------
% SpikeInterface wants data in 'double' format
%------------------------------------------------------------------------
exportOpts.OutputFormat = 'double';

%{
from https://spikeinterface.readthedocs.io/en/latest/how_to/load_matlab_data.html

----------------------------------------------------
Writing data for SpikeInterface:
----------------------------------------------------
% Define the size of your data
numSamples = 1000;
numChannels = 384;

% Generate random data as an example
data = rand(numSamples, numChannels);

% Write the data to a binary file
fileID = fopen('your_data_as_a_binary.bin', 'wb');
fwrite(fileID, data, 'double');
fclose(fileID);
%}


%------------------------------------------------------------------------
%% run!
%------------------------------------------------------------------------
[nD, nI, nW] = export_raw(exportOpts);


%------------------------------------------------------------------------
%% as a check, read in binary data
%------------------------------------------------------------------------
nC = length(exportOpts.Channels);

fp = fopen(fullfile(exportOpts.OutputPath, exportOpts.OutputFile), ...
                  'rb');
M = fread(fp, exportOpts.OutputFormat);
fclose(fp);

sM = size(M);
fprintf('read %d elements\n', numel(M));
fprintf('Size of data read from file: [%d, %d]\n', sM(1), sM(2));
M2 = reshape(M, [numel(M)/nC, nC]);
size(M2)


%{
%% do a simple test
nC = length(exportOpts.Channels);

a = ones(nC, 20);
for n = 1:nC
   a(n, :) = n*a(n, :);
end

%%

% open raw file
fp = fopen('test.bin', 'wb');
nw = fwrite(fp, a, exportOpts.OutputFormat);
fclose(fp);

% read raw file
fp = fopen('test.bin', 'rb');
b  = fread(fp, exportOpts.OutputFormat);
fclose(fp);
%%
reshape(b, size(a))

reshape(b, [nC, numel(b)/nC])
%}

