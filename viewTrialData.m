function varargout = viewTrialData(varargin)
%------------------------------------------------------------------------
% viewTrialData(datfile)
%------------------------------------------------------------------------
% optosort project
%------------------------------------------------------------------------
% reads in data from .dat files created by the opto program
% Then opens common_ref_viewer() to allow user to page through data trials
% and view 16 channels of recorded data. 
% User can apply common average reference or common median reference to
% trial data to explore effects on recording
%
% Limitations: 
%  .dat files must have 16 channels of data recorded
%------------------------------------------------------------------------
% Input Args:
%  datfile     (optional) filename for .dat file (with path prepended)
%
% Output Args: 
%  figure handle
%------------------------------------------------------------------------
% See also: common_ref_viewer, opto program
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Created: 2 May, 2023 Sharad Shanbhag  sshanbhag@neomed.edu
% 
% Revisions:
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% settings
%------------------------------------------------------------------------
% read in all 16 channels
Channels = 1:16;

%------------------------------------------------------------------------
% filter parameters for raw neural data
%------------------------------------------------------------------------
% [highpass lowpass] cutoff frequencies in Hz
BPfilt.Fc = [250 4000];
% order of filter. note that the filtfilt() function in MATLAB is used,
% so the effective order is doubled. typically use 5
BPfilt.forder = 5;
% ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% transients from filtering
BPfilt.ramp = 1;
% filter type. 'bessel' or 'butter' (for butterworth). typically use butter
% as bessel is low-pass only and we need to remove low frequency crap from
% the raw data
BPfilt.type = 'butter';
% don't resample data (pass in empty value)
resampleData = [];

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Data file - provided as input or get from user
%------------------------------------------------------------------------
if isempty(varargin)
   % get .dat file from user 
   [datfile, datpath] = uigetfile('*.dat', ...
                                    'Open .dat file from opto program');

   % if fname is 0 then user hit cancel button, exit function
   if isequal(datfile, 0)
      return
   end
else
   % make sure input file exists
   if ~exist(varargin{1}, 'file')
      error('%s: datafile not found %s', mfilename, varargin{1});
   end
   [datpath, datfile, dext] = fileparts(varargin{1});
   datfile = [datfile dext];
end

%------------------------------------------------------------------------
% build _testdata.mat filename
%------------------------------------------------------------------------
[~, tbase] = fileparts(datfile);
tname = [tbase '_testdata.mat'];
% check that it exists
if ~exist(fullfile(datpath, tname), 'file')
   error('%s: _testdata.mat file not found %s', fullfile(datpath, tname));
end
% define path to data file
F = defineSampleData({datpath}, {datfile}, {tname});

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% read and filter data, no resampling
%------------------------------------------------------------------------
% get the data and information from the raw files
%	cSweeps is a {nfiles, 1} cell array with each element holding an
% 	{nChannels X ntrials} cell array of data for each file
%	nexInfo is an instance of a SpikeInfo object
%%------------------------------------------------------------------------
[cSweeps, nexInfo] = read_data_for_export(F, Channels, ...
                                                BPfilt, resampleData);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% call dataExplore with sweep data and nexInfo
%------------------------------------------------------------------------
%------------------------------------------------------------------------
cfig = common_ref_viewer(cSweeps, nexInfo);

varargout{1} = cfig;

