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
F.DataPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort';
% F.DataPath = '/Volumes/SJS_XFER/Work/MT-IC-R-data';
%------------------------------------------------------------------------
% DATA FILES
%------------------------------------------------------------------------
% F.DataFile = {'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
% 				'1372_20191126_03_01_1500_BBN.dat'; ...
% 				'1372_20191126_03_01_1500_FRA.dat'; ...
% 				'1372_20191126_03_01_1500_WAV.dat'; };
% F.TestFile = {'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
% 				'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
% 				'1372_20191126_03_01_1500_FRA_testdata.mat'; ...
% 				''; };
F.DataFile = {'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
				'1372_20191126_03_01_1500_BBN.dat'; ...
				'1372_20191126_03_01_1500_FRA.dat'; };
F.TestFile = {'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
				'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
				'1372_20191126_03_01_1500_FRA_testdata.mat'; };
F.Channels = [11 9 14];

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
% run!
%------------------------------------------------------------------------
nD = export_for_plexon(F);

