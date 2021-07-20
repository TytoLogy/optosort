%------------------------------------------------------------------------
% export_working_infusion.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% sample script - working file for Sharad during development of infusion
% data analysis
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 19 July. 2021 (SJS)
%
% Revisions:
%	
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
% 1459 Opto data
%---------------------------------------
exportOpts.testData = false;
% this is for data on SJS's Linux machine
exportOpts.DataPath = '/media/Data/NeuroData/Raw/1466/20210617';
% local working dir
% exportOpts.OutputPath = '/home/sshanbhag/Work/Data/Test';
% on myotis (ubuntu linux desktop):
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/sort/1466';
% on NAS:
% exportOpts.OutputPath = '/media/NAS/By User/SJS/tmp'
exportOpts.DataFile = {	...
   '1466_20210617_01_0_550_BBN.dat', ...
   '1466_20210617_01_0_550_FREQ_TUNING.dat', ...
   '1466_20210617_01_0_550_FRA.dat', ...
   '1466_20210617_01_0_550_WAV.dat', ...
   '1466_20210617_01_0_550_CLICK.dat', ...
   '1466_20210617_01_0_550_BBN_5.dat', ...
   '1466_20210617_01_0_550_BBN_10.dat', ...
   '1466_20210617_01_0_550_BBN_15.dat', ...
   '1466_20210617_01_0_550_BBN_20.dat', ...
   '1466_20210617_01_0_550_BBN_postAnest.dat', ...
   '1466_20210617_01_0_550_FREQ_TUNING_postAnest.dat', ...
   '1466_20210617_01_0_550_FRA_postAnest.dat', ...
   '1466_20210617_01_0_550_WAV_postAnest.dat', ...
   '1466_20210617_01_0_550_CLICK_postAnest.dat', ...
   };
exportOpts.TestFile = { ...
   '1466_20210617_01_0_550_BBN_testdata.mat', ...
   '1466_20210617_01_0_550_FREQ_TUNING_testdata.mat', ...
   '1466_20210617_01_0_550_FRA_testdata.mat', ...
   '1466_20210617_01_0_550_WAV_testdata.mat', ...
   '1466_20210617_01_0_550_CLICK_testdata.mat', ...
   '1466_20210617_01_0_550_BBN_5_testdata.mat', ...
   '1466_20210617_01_0_550_BBN_10_testdata.mat', ...
   '1466_20210617_01_0_550_BBN_15_testdata.mat', ...
   '1466_20210617_01_0_550_BBN_20_testdata.mat', ...
   '1466_20210617_01_0_550_BBN_postAnest_testdata.mat', ...
   '1466_20210617_01_0_550_FREQ_TUNING_postAnest_testdata.mat', ...
   '1466_20210617_01_0_550_FRA_postAnest_testdata.mat', ...
   '1466_20210617_01_0_550_WAV_postAnest_testdata.mat', ...
   '1466_20210617_01_0_550_CLICK_postAnest_testdata.mat', ...
   };
exportOpts.OutputFile = '1466_20210617_01_0_550_MERGEbbn_C1-8.nex';
% specify which events to write
% Valid event names:
%    'all'             write all events
%    'startsweep'      start of each sweep (aka trial, rep)
%    'endsweep'        end of each sweep
%    'filestart'       beginning of each file's data
%    'fileend'         end of each file's data
%    'filename'        event with source .dat file name at start of 
%                      that files' data
%    'stimstart'       stimulus onset
%    'stimend'         stimulus end
%    'stim_id'         name & properties of stimulus at onset
exportOpts.eventsToWrite = {'filename', 'startsweep', 'stimstart'};
exportOpts.Channels = 1:8;

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
