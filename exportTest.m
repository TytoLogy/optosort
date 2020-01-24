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
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

clear all

if ~exist('readOptoData', 'file')
	addOptoPaths
end

%------------------------------------------------------------------------
% PATHS TO DATA FILES
%------------------------------------------------------------------------
% can specify individual file paths... these must match the paths in
% F.DataFile!!!!!!
% F.DataPath = {'~/Work/Data/TestData/MT_IC';
% 				'~/Work/Data/TestData/MT_IC'; ...
% 				'~/Work/Data/TestData/MT_IC'; };
% or single path
% F.DataPath = '~/Work/Data/TestData/MT_IC';
% F.DataPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort';
% F.DataPath = '/Volumes/SJS_XFER/Work/MT-IC-R-data';
%------------------------------------------------------------------------
% DATA FILES
%------------------------------------------------------------------------

%---------------------------------------
% IC data from probe, include WAV test
% F.DataPath = ['/Volumes/Wenstrup Laboratory/By User/SJS/Data' ...
% 							'/SpikeSort/IC-probe/1372'];
% F.DataFile = {'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
% 				'1372_20191126_03_01_1500_BBN.dat'; ...
% 				'1372_20191126_03_01_1500_FRA.dat'; ...
% 				'1372_20191126_03_01_1500_WAV.dat'; };
% F.TestFile = {'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
% 				'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
% 				'1372_20191126_03_01_1500_FRA_testdata.mat'; ...
% 				''; };
% F.Channels = [11 9 14];
%---------------------------------------
% IC data from probe, NO WAV test
% F.DataPath = ['/Volumes/Wenstrup Laboratory/By User/SJS/Data' ...
% 							'/SpikeSort/IC-probe/1372'];
F.DataPath = '~/Work/Data/TestData/MT-IC-R-data';
F.DataFile = {'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
				'1372_20191126_03_01_1500_BBN.dat'; ...
				'1372_20191126_03_01_1500_FRA.dat'; };
F.TestFile = {'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
				'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
				'1372_20191126_03_01_1500_FRA_testdata.mat'; };
F.Channels = [11 9 14];
F.OutputPath = F.DataPath;
F.OutputFile = '1372_20191126_03_01_1500_test.nex';
%---------------------------------------

%---------------------------------------
% MG sample data with probe
% 1382_20191212_02_02_3200_WAV.dat
% 1382_20191212_02_02_3200_WAV_PSTH.fig
% 1382_20191212_02_02_3200_WAV_wavinfo.mat
% F.DataPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort/MG';
% F.DataFile = {	'1382_20191212_02_02_3200_FREQ_TUNING.dat'; ...
% 					'1382_20191212_02_02_3200_BBN.dat'; ...
% 					'1382_20191212_02_02_3200_FRA.dat';	};
% F.TestFile = {	'1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'; ...
% 					'1382_20191212_02_02_3200_BBN_testdata.mat'; ...
% 					'1382_20191212_02_02_3200_FRA_testdata.mat'; };
% F.Channels = [4, 5, 11, 14];
%---------------------------------------

%---------------------------------------
% IC sample data tungsten (single channel)
% 1323_20190722_03_02_632_WAV.dat
% 1323_20190722_03_02_632_WAV_testdata.mat
% 1323_20190722_03_02_632_WAV_wavinfo.mat
% F.DataPath = ['/Volumes/Wenstrup Laboratory/By User/SJS/Data' ...
% 						'/SpikeSort/IC-tungsten/'];
% F.DataFile = {	'1323_20190722_03_02_632_FREQ_TUNING.dat'; ...
% 					'1323_20190722_03_02_632_BBN.dat'; ...
% 					'1323_20190722_03_02_632_TONE_LEVEL.dat';	};
% F.TestFile = {	'1323_20190722_03_02_632_FREQ_TUNING_testdata.mat'; ...
% 					'1323_20190722_03_02_632_BBN_testdata.mat'; ...
% 					'1323_20190722_03_02_632_TONE_LEVEL_testdata.mat'; };
% F.Channels = [8];



%------------------------------------------------------------------------
% filter parameters for raw neural data
%------------------------------------------------------------------------
% [highpass lowpass] cutoff frequencies in Hz
BPfilt.Fc = [300 4000];
% order of filter. note that the filtfilt() function in MATLAB is used,
% so the effective order is doubled. 
BPfilt.forder = 5;
% ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% transients from filtering
BPfilt.ramp = 1;
% assign filter spec to F struct
F.BPfilt = BPfilt;

%------------------------------------------------------------------------
%% run!
%------------------------------------------------------------------------
nD = export_for_plexon(F);

