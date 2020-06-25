%------------------------------------------------------------------------
% plot_curve_demo.m
%------------------------------------------------------------------------
% TytoLogy:optosort
%--------------------------------------------------------------------------
% example script for generating plots of data
%------------------------------------------------------------------------
% See Also: optoproc, opto (TytoLogy:opto program)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 25 June 2020 (SJS)
%	 
% Revisions:
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%% add path to readPLXFileC if needed
%------------------------------------------------------------------------
% readPLXFileC is a function downloaded from MATLAB Central that allows
% direct reading of PLX file data in Matlab
if ~exist('readPLXFileC', 'file')
	fprintf('plot_curve_demo: adding readPLXFile to path\n');
	addpath('readPLXFileC');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% define paths and data files
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% data locations
plxFilePath = '/Users/sshanbhag/Work/Data/TestData/working';
nexInfoPath = plxFilePath;

% sorted data file
plxFile = '1407_20200309_03_01_1350_MERGEALL.plx';

% nexinfo file
nexInfoFile = '1407_20200309_03_01_1350_MERGEALL_nexinfo.mat';

sendmsg(sprintf('Using data from file: %s', plxFile));
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% load sorted data
%------------------------------------------------------------------------
% How to use:
% import_from_plexon(<plx file name>, <nexinfo file name>, 
%								<'continuous'/'nocontinuous'>)
%
%	will return a SpikeData object containing data from plx file:
% 	- spike times/unit information if sorted
% 	- continuous data, if saved in plx and 'continuous' is specified (default)
% 	- file stimulus info
%------------------------------------------------------------------------
%------------------------------------------------------------------------
S = import_from_plexon(fullfile(plxFilePath, plxFile), ...
							fullfile(nexInfoPath, nexInfoFile), 'continuous');
						
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% information about file
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% show file, channel, unit
%------------------------------------------------------------------------
% display file, channel, unit info
[fileList, channelList, unitList] = S.printInfo;

%------------------------------------------------------------------------
% Show how to list files and curve types
%------------------------------------------------------------------------
% loop through the number of files merged into file to be sorted
sendmsg('File and test information:')
for f = 1:S.Info.nFiles
	% for each file number, display the test type and test name (i know why
	% two different things? an effect of an old kludge...)
	fprintf('File %d:\n', f);
	fprintf('\tTest Type: %s\n', S.Info.FileInfo{f}.testtype);
	fprintf('\tTest Name: %s\n', S.Info.FileInfo{f}.testname);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot freq-tuning curves for non-zero units
%------------------------------------------------------------------------
%------------------------------------------------------------------------






