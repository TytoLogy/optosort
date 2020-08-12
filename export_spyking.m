%------------------------------------------------------------------------
% export_spyking.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% exports data in raw format for use with spyKing CIRCUS sorting program
% sample script - working file for Sharad during development
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 13 July, 2020 (SJS)
%	from export_working
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
% 1429 Cambridge NeuroTech Probe
%---------------------------------------
exportOpts.DataPath = '/Users/sshanbhag/Work/Data/TestData/working/1429/20200707';
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/working/exports/1429/spyk';
exportOpts.DataFile = '1429_20200707_01_01_2942_BBN.dat';
exportOpts.OutputFile = '1429_20200707_01_01_2942_BBN.bin';
exportOpts.TestFile = '1429_20200707_01_01_2942_BBN_testdata.mat';
exportOpts.Channels = 1:16;
exportOpts.testData = false;
%--------------------------------------
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

if ~exist(exportOpts.OutputPath)
	mkdir(exportOpts.OutputPath);
end
%------------------------------------------------------------------------
% resample data to nearest lower integer value?
%------------------------------------------------------------------------
% exportOpts.resampleData = [];

%------------------------------------------------------------------------
% run!
%------------------------------------------------------------------------
[nD, nI, nP] = export_raw(exportOpts);