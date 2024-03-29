%------------------------------------------------------------------------
% export_working.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% sample script - working file for Sharad during development
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 8 January, 2020 (SJS)
%
% Revisions:
%	26 Mar 2020 (SJS): reworking for use
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% check to make sure tytoLogy (esp opto) things are on path
if ~exist('readOptoData', 'file')
	fprintf(['readOptoData (and possibly other TytoLogy/opto functions' ...
						' not found!\n'])
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
% exportOpts.DataPath = '~/Work/Data/TestData/MT/1382';
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
% 1407
%---------------------------------------
%{
exportOpts.DataPath = '~/Work/Data/TestData/MT/1407';
exportOpts.DataFile = {	'1407_20200309_03_01_1350_FREQ_TUNING.dat'; ...
					'1407_20200309_03_01_1350_BBN.dat'; ...
					'1407_20200309_03_01_1350_FRA.dat';	};
exportOpts.TestFile = {	'1407_20200309_03_01_1350_FREQ_TUNING_testdata.mat'; ...
					'1407_20200309_03_01_1350_BBN_testdata.mat'; ...
					'1407_20200309_03_01_1350_FRA_testdata.mat'; };
exportOpts.Channels = [4, 5, 7, 15];
exportOpts.OutputPath = exportOpts.DataPath;
exportOpts.OutputFile = '1407_20200309_03_01_1350_MERGEEVENTS.nex';
%}
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
%---------------------------------------
% 1407 multichannel, all files
%---------------------------------------

exportOpts.DataPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
% write to D drive on PETROL for testing with Plexon OFS
% exportOpts.OutputPath = '/Volumes/D/1407';
% local working dir
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/working';
exportOpts.DataFile = {	'1407_20200309_03_01_1350_BBN.dat', ...
								'1407_20200309_03_01_1350_FREQ_TUNING.dat', ...
								'1407_20200309_03_01_1350_WAV.dat', ...
								'1407_20200309_03_01_1350_FRA.dat'	};
exportOpts.TestFile = { '1407_20200309_03_01_1350_BBN_testdata.mat', ...
								'1407_20200309_03_01_1350_FREQ_TUNING_testdata.mat', ...
								'1407_20200309_03_01_1350_WAV_testdata.mat', ...
								'1407_20200309_03_01_1350_FRA_testdata.mat'	};
exportOpts.OutputFile = '1407_20200309_03_01_1350_MERGEEVENTS.nex';
exportOpts.Channels = [4, 5, 7, 15];
exportOpts.testData = false;
%--------------------------------------
%---------------------------------------
%}

%{
%---------------------------------------
% 1407 multichannel, test with smaller files
%---------------------------------------
exportOpts.DataPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
% write to D drive on PETROL for testing with Plexon OFS
% exportOpts.OutputPath = '/Volumes/D/1407';
% local working dir
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/working';
exportOpts.DataFile = {	'1407_20200309_03_01_1350_BBN.dat', ...
								'1407_20200309_03_01_1350_FREQ_TUNING.dat' 	};
exportOpts.TestFile = { '1407_20200309_03_01_1350_BBN_testdata.mat', ...
								'1407_20200309_03_01_1350_FREQ_TUNING_testdata.mat'	};
exportOpts.OutputFile = '1407_20200309_03_01_1350_MERGEEVENTS_BBNFREQ.nex';
exportOpts.Channels = [4, 5, 7, 15];
exportOpts.testData = false;
%---------------------------------------
%---------------------------------------
%}

%{
%---------------------------------------
% 1407 FRA Only
%---------------------------------------
exportOpts.DataPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
% write to D drive on PETROL for testing with Plexon OFS
% exportOpts.OutputPath = '/Volumes/D/1407';
% local working dir
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/working';
exportOpts.DataFile = {	'1407_20200309_03_01_1350_FRA.dat'	};
exportOpts.TestFile = { '1407_20200309_03_01_1350_FRA_testdata.mat'	};
exportOpts.OutputFile = '1407_20200309_03_01_1350_FRA.nex';
exportOpts.Channels = [4, 5, 7, 15];
exportOpts.testData = false;
%}

%{
%---------------------------------------
% 1458 Opto data
%---------------------------------------
% this is for data on SJS's Linux machine
exportOpts.DataPath = '/media/Data/NeuroData/Raw/1458/20210506';
% local working dir
% exportOpts.OutputPath = '/home/sshanbhag/Work/Data/Test';
% on myotis (ubuntu linux desktop):
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/sort/1458';
% on NAS:
% exportOpts.OutputPath = '/media/NAS/By User/SJS/tmp'
%{
exportOpts.DataFile = {	'1458_20210506_01_0_3300_BBN.dat', ...
                        '1458_20210506_01_0_3300_FREQ_TUNING.dat', ...
                        '1458_20210506_01_0_3300_FRA.dat', ...
                        '1458_20210506_01_0_3300_WAV.dat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-5ms.dat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-30ms.dat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-100ms.dat', ...
                        ...
                        };
exportOpts.TestFile = { '1458_20210506_01_0_3300_BBN_testdata.mat', ...
                        '1458_20210506_01_0_3300_FREQ_TUNING_testdata.mat', ...
                        '1458_20210506_01_0_3300_FRA_testdata.mat', ...
                        '1458_20210506_01_0_3300_WAV_testdata.mat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-5ms_testdata.mat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-30ms_testdata.mat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-100ms_testdata.mat', ...
                        ...
                        };
exportOpts.OutputFile = '1458_20210506_01_0_3300_MERGE_3opt.nex';
%}
exportOpts.DataFile = {	'1458_20210506_01_0_3300_BBN.dat', ...
                        '1458_20210506_01_0_3300_FREQ_TUNING.dat', ...
                        '1458_20210506_01_0_3300_FRA.dat', ...
                        '1458_20210506_01_0_3300_WAV.dat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-100ms.dat', ...
                        '1458_20210506_01_0_3300_CLICK.dat', ...
                        ...
                        };
exportOpts.TestFile = { '1458_20210506_01_0_3300_BBN_testdata.mat', ...
                        '1458_20210506_01_0_3300_FREQ_TUNING_testdata.mat', ...
                        '1458_20210506_01_0_3300_FRA_testdata.mat', ...
                        '1458_20210506_01_0_3300_WAV_testdata.mat', ...
                        '1458_20210506_01_0_3300_OPTO-AMP-100ms_testdata.mat', ...
                        '1458_20210506_01_0_3300_CLICK_testdata.mat', ...
                        ...
                        };
exportOpts.OutputFile = '1458_20210506_01_0_3300_MERGE_sgloptoclick.nex';
exportOpts.Channels = [8, 10, 15];
exportOpts.testData = false;
%}

%---------------------------------------
% 1459 Opto data
%---------------------------------------
exportOpts.testData = false;
% this is for data on SJS's Linux machine
exportOpts.DataPath = '/media/Data/NeuroData/Raw/1459/20210422';
% local working dir
% exportOpts.OutputPath = '/home/sshanbhag/Work/Data/Test';
% on myotis (ubuntu linux desktop):
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/sort/1459';
% on NAS:
% exportOpts.OutputPath = '/media/NAS/By User/SJS/tmp'
exportOpts.DataFile = {	'1459_20210422_02_0_3500_BBN.dat', ...
                        '1459_20210422_02_0_3500_FREQ_TUNING.dat', ...
                        '1459_20210422_02_0_3500_FRA.dat', ...
                        '1459_20210422_02_0_3500_WAV.dat', ...
                        '1459_20210422_02_0_3500_OPTO-AMP.dat', ...
                        '1459_20210422_02_0_3500_CLICK.dat', ...
                        ...
                        };
exportOpts.TestFile = { '1459_20210422_02_0_3500_BBN_testdata.mat', ...
                        '1459_20210422_02_0_3500_FREQ_TUNING_testdata.mat', ...
                        '1459_20210422_02_0_3500_FRA_testdata.mat', ...
                        '1459_20210422_02_0_3500_WAV_testdata.mat', ...
                        '1459_20210422_02_0_3500_OPTO-AMP_testdata.mat', ...
                        '1459_20210422_02_0_3500_CLICK_testdata.mat', ...
                        ...
                        };
exportOpts.OutputFile = '1459_20210422_02_0_3500_MERGE_allchan.nex';
exportOpts.Channels = 1:16;
% exportOpts.Channels = [6, 7, 14, 16];
%{
% just export opto for now to see what channels to export for real
exportOpts.DataFile = {	'1459_20210422_02_0_3500_OPTO-AMP.dat' ...
                        };
exportOpts.TestFile = { '1459_20210422_02_0_3500_OPTO-AMP_testdata.mat' ...
                        ...
                        };
exportOpts.OutputFile = '1459_20210422_02_0_3500_MERGE_optoonly.nex';
exportOpts.Channels = 1:16;
%}

%---------------------------------------
%---------------------------------------
%------------------------------------------------------------------------
% filter parameters for raw neural data
%------------------------------------------------------------------------
% [highpass lowpass] cutoff frequencies in Hz
exportOpts.BPfilt.Fc = [250 4000];
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

%------------------------------------------------------------------------
% resample data to nearest lower integer value?
%------------------------------------------------------------------------
% exportOpts.resampleData = [];

%------------------------------------------------------------------------
% run!
%------------------------------------------------------------------------
[nD, nI] = export_for_plexon(exportOpts);

% save('testobj.mat', 'nD', 'nI', '-MAT')
