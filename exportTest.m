%------------------------------------------------------------------------
% exportTest.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% sample script to show how to export data for spike sorting
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
						'not found!\n'])
	fprintf('Please check Matlab path(s)\n');
	fprintf(['e.g.:\n' ...
				'addpath(''~/Work/Code/Matlab/dev/TytoLogy/Experiments/Opto'')\n']);
end
	
% add the opto objects path
addpath('OptoObjects');

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
exportOpts.DataPath = '~/Work/Data/TestData/MT';
exportOpts.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'; ...
					'1382_20191212_02_02_3200_BBN.dat'; ...
					'1382_20191212_02_02_3200_FRA.dat';	};
exportOpts.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'; ...
					'1382_20191212_02_02_3200_BBN_testdata.mat'; ...
					'1382_20191212_02_02_3200_FRA_testdata.mat'; };
exportOpts.Channels = [4, 5, 11, 14];
exportOpts.OutputPath = exportOpts.DataPath;
exportOpts.OutputFile = '1382_20191212_02_02_3200_MERGE.nex';
% exportOpts.OutputFile = '1382_20191212_02_02_3200_CURVEOBJ.nex';
% exportOpts.OutputPath = pwd;
% exportOpts.OutputFile = 'test.nex';
%---------------------------------------


%------------------------------------------------------------------------
% filter parameters for raw neural data
%------------------------------------------------------------------------
% [highpass lowpass] cutoff frequencies in Hz
exportOpts.BPfilt.Fc = [300 4000];
% order of filter. note that the filtfilt() function in MATLAB is used,
% so the effective order is doubled. typically use 5
exportOpts.BPfilt.forder = 5;
% ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% transients from filtering
exportOpts.BPfilt.ramp = 1;
% filter type. 'bessel' or 'butter' (for butterworth). typically use bessel
exportOpts.BPfilt.type = 'bessel';

%------------------------------------------------------------------------
% run!
%------------------------------------------------------------------------
[nD, nI] = export_for_plexon(exportOpts);

