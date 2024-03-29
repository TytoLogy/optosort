%------------------------------------------------------------------------
% click_working.m
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
% Created: 24 June, 2021 (SJS)
%  code to develop export/processing of click stimulus data
% Revisions:
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

%{
%---------------------------------------
%---------------------------------------
% 1458 Opto, Click data
%---------------------------------------
% use this for working with the older click data files without stimList
%---------------------------------------
% this is for data on SJS's Linux machine
exportOpts.DataPath = '/media/Data/NeuroData/Raw/1458/20210506';

% write to D drive on PETROL for testing with Plexon OFS
% exportOpts.OutputPath = '/Volumes/D/1407';
% local working dir
% exportOpts.OutputPath = '/home/sshanbhag/Work/Data/Test';
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/click';
exportOpts.DataFile = {	%'1458_20210506_01_0_3300_BBN.dat', ...
								%'1458_20210506_01_0_3300_OPTO-AMP-100ms.dat', ...
                        '1458_20210506_01_0_3300_CLICK.dat', ...
                        ...
                        };
exportOpts.TestFile = { ...
              %'1458_20210506_01_0_3300_BBN_testdata.mat', ...
              %'1458_20210506_01_0_3300_OPTO-AMP-100ms_testdata.mat', ...
              '1458_20210506_01_0_3300_CLICK_testdata.mat', ...
                  ...
                  };
exportOpts.OutputFile = '1458_20210506_01_0_3300_oldclick_MERGE.nex';
exportOpts.Channels = [8];
exportOpts.testData = false;
%}


%---------------------------------------
%---------------------------------------
% 000 Opto, Click data
%---------------------------------------
% use this for working with the newer 
% click data files, at least until real
% click data in the new format are available
%---------------------------------------
% this is for data on SJS's Linux machine
exportOpts.DataPath = '/media/Data/NeuroData/TestData/click/20210623';

% write to D drive on PETROL for testing with Plexon OFS
% exportOpts.OutputPath = '/Volumes/D/1407';
% local working dir
% exportOpts.OutputPath = '/home/sshanbhag/Work/Data/Test';
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/click';
exportOpts.DataFile = {	%'1458_20210506_01_0_3300_BBN.dat', ...
								%'1458_20210506_01_0_3300_OPTO-AMP-100ms.dat', ...
                        '000_20210623_0_0_0_CLICK.dat', ...
                        ...
                        };
exportOpts.TestFile = { ...
              %'1458_20210506_01_0_3300_BBN_testdata.mat', ...
              %'1458_20210506_01_0_3300_OPTO-AMP-100ms_testdata.mat', ...
              '000_20210623_0_0_0_CLICK.dat_testdata.mat', ...
                  ...
                  };
exportOpts.OutputFile = '000_20210623_0_0_0_newclick_MERGE.nex';
exportOpts.Channels = [8];
exportOpts.testData = false;
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
exportOpts.resampleData = [];

%------------------------------------------------------------------------
% run export_for_plexon!
%------------------------------------------------------------------------
[nD, nI] = export_for_plexon(exportOpts);

ci = nI.FileInfo{1}
% 
% save('testobj.mat', 'nD', 'nI', '-MAT')



%%

if isempty(ci.Dinf.stimList)
   % no stimList, so we need to make some assumptions
   % get varied values
   levels = ci.varied_values;
   % check if NULL stimulus was played
   if ci.Dinf.test.NullStim
      % since null stim was played, ASSUME that stimindex == 1 is null,
      % stimindex == 2 is click
      % create list of levels by sweep
      levellist = zeros(size(ci.Dinf.test.stimIndices));
      % set levels for stimIndices == 1(null) to 0;
      levellist(ci.Dinf.test.stimIndices == 1) = levels(1);
      % set levels for stimIndices == 2 (click) to level
      levellist(ci.Dinf.test.stimIndices == 2) = levels(2);      
      
   end
end
