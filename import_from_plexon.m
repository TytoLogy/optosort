% function varargout = import_from_plexon(varargin)
%------------------------------------------------------------------------
% import_from_plexon.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% working script for importing sorted data from plexon
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
% Initial things to define
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
sortedPath = '~/Work/Data/TestData/MT_IC';
rawPath = '~/Work/Data/TestData/MT_IC';
nexPath = '~/Work/Data/TestData/MT_IC';

% sorted data file
% sortedFile = '1323_20190722_03_02_632_MERGE.mat';
sortedFile = '1382_20191212_02_02_3200.mat';

% nexinfo file
% nexInfoFile = '1323_20190722_03_02_632_MERGE_nexinfo.mat';
nexInfoFile = '1382_20191212_02_02_3200_MERGE_nexinfo.mat';

% nex file
nexFile = '1382_20191212_02_02_3200_MERGE.mat';

%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n', sepstr);
fprintf('import_from_plexon running...\n');
fprintf('%s\n', sepstr);


%------------------------------------------------------------------------
%% load data
%------------------------------------------------------------------------

% nexinfo
load(fullfile(nexPath, nexInfoFile), 'nexInfo');

% plexon sorted data
plxvars = who('-file', fullfile(sortedPath, sortedFile));
if isempty(plxvars)
	error('No variables in plexon output file %s', ...
					fullfile(sortedPath, sortedFile));
else
	fprintf('\n%s\n', sepstr);
	fprintf('Found %d channels in %s\n', length(plxvars), ...
							fullfile(sortedPath, sortedFile))
	fprintf('\n%s\n', sepstr);
	for n = 1:length(plxvars)
		fprintf('%s\n', plxvars{n});
	end
	fprintf('%s\n', sepstr);
end	
load(fullfile(sortedPath, sortedFile))

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Plexon-exported sorted data:
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Column 1: unit number (where 0 is unsorted)
% Column 2: timestamp where spike crosses threshold (in seconds)
% Columns 3-34 (assuming waveform window of 1311us / 32 samples):
%   waveform snippet, with or without prewindow as set in Offline Sorter
%   (prewindow default: 494us / 12 samples)
%   (window default: 1311us / 32 samples)
%   This is in units of samples/sec of raw data file 
%   (24414.063 Hz based on settings in data acquisition program 
%        HPSearch or PresentStimCurve in RosenLab)

% may not be accurate.... might be
% Column 1: channel (?)
% Column 2: unit #
% Column 3: timestamp (in seconds)
% Column 4-... : waveform
%
% 22 Jan 2019
% if waveforms are clipped, look at pp 109 and 118 in offline sorter manual
% to fix gain
% 
% other issue: transients on first 1-4 samples of each waveform.... ????
% just noise????
%------------------------------------------------------------------------

