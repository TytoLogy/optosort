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

% sorted data file
sortedFile = '1372_20191126_03_01_1500_FREQ_TUNING_plx.mat';

% nexinfo file
nexInfoFile = '1372_20191126_03_01_1500_nexinfo.mat';

%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n', sepstr);
fprintf('import_from_plexon running...\n');
fprintf('%s\n', sepstr);


%------------------------------------------------------------------------
% load data
%------------------------------------------------------------------------
% plexon sorted data

% nexinfo
load(nexInfoFile)

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

