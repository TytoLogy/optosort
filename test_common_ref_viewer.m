%------------------------------------------------------------------------
% test_common_ref_viewer
%------------------------------------------------------------------------
% Script to load .dat file and then call common_ref_viewer to plot
% recording and test effects of common average and median referencing
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% load data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% specify file to examine
%  path, dat file, _testdata.mat file
%------------------------------------------------------------------------
if strcmpi(computer, 'GLNXA64')
   dpath = '/media/sshanbhag/SSData/Data/Mouse/Raw/BLA/1500/20230213';
else
   dpath = '/Volumes/SSData/Data/Mouse/Raw/BLA/1500/20230213';
end
% dname = '1500_20230213_01_0_3352_BBN.dat';
% tname = '1500_20230213_01_0_3352_BBN_testdata.mat';
dname = '1500_20230213_01_0_3352_WAV.dat';
tname = '1500_20230213_01_0_3352_WAV_testdata.mat';

% define path to data file
F = defineSampleData({dpath}, {dname}, {tname});

%------------------------------------------------------------------------
% get the data and information from the raw files
%	cSweeps is a {nfiles, 1} cell array with each element holding an
% 	{nChannels X ntrials} cell array of data for each file
%	nexInfo is an instance of a SpikeInfo object
%
% use empty values for BPfilt (arg 3) and resampleData (arg 4) so 
% that no filtering or resampling is done
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
%% read and filter data, no resampling
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[cSweeps, nexInfo] = read_data_for_export(F, Channels, ...
                                                BPfilt, resampleData);
                

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% call dataExplore with sweep data and nexInfo
%------------------------------------------------------------------------
%------------------------------------------------------------------------
cfig = common_ref_viewer(cSweeps, nexInfo);



