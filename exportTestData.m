%------------------------------------------------------------------------
% exportTestData.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% generates 
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 9 June, 2020 (SJS)
%
% Revisions:
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



%---------------------------------------
%---------------------------------------
% testing scenarios
%---------------------------------------


% 1407 multichannel
exportOpts.DataPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/working';
exportOpts.DataFile = '1407_20200309_03_01_1350_BBN.dat';
exportOpts.OutputFile = '1407_20200309_03_01_1350_TIMETESTDATA.nex';
exportOpts.TestFile = '1407_20200309_03_01_1350_BBN_testdata.mat';
exportOpts.Channels = 1;

%---------------------------------------
%---------------------------------------

%{
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
% filter type. 'bessel' or 'butter' (for butterworth). typically use butter
% as bessel is low-pass only and we need to remove low frequency crap from
% the raw data
exportOpts.BPfilt.type = 'butter';
%}

%------------------------------------------------------------------------
% resample data to nearest lower integer value?
%------------------------------------------------------------------------
exportOpts.resampleData = [];

%------------------------------------------------------------------------
% run!
%------------------------------------------------------------------------
[nD, nI] = export_for_plexon(exportOpts);

