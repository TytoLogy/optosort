%------------------------------------------------------------------------
% optosort_working.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% working script for opto data <--> plexon
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

%% Initial things

sepstr = '----------------------------------------------------';
fprintf('\n%s\n', sepstr);
fprintf('%s running...\n', mfilename);
fprintf('%s\n', sepstr);

% add path to .nex file utils
if ~exist('nexCreateFileData', 'file')
	fprintf('Adding NEX file utility paths\n')
	addpath(['~/Work/Code/Matlab/stable/Toolbox/NeuroExplorer' ...
					'/HowToReadAndWriteNexAndNex5FilesInMatlab']);
end
% add path to opto if needed
if ~exist('readOptoData', 'file')
	fprintf('Adding opto paths\n')
	addOptoPaths
end

% define path to data file and data file for testing
defineSampleData
fprintf('DataPath = %s\n', DataPath);
fprintf('DataFile = %s\n', DataFile);
fprintf('Animal: %s\n', F.animal);

% !!!!!!!
% in future, will need to have user specify all the files for this
% recording session and recording location!!!
% ¡¡¡¡¡¡¡

% channel(s) of data to obtain
% channel = 8;
channel = [1 2 3];

%% read in data

% can probably just use readOptoData
% % the data will be further processed during sorting. for now, use a fairly
% % broad filter
% filtband = [5000 10000];
% % get the data - only need the raw data (stored in D cell array of each
% % stimulus presentation) and the information in Dinf about the test.
% [D, Dinf] = getFilteredOptoData(fullfile(datadir, datafile), ...
% 											'filter', filtband, ...
% 											'channel', channel);
										
[D, Dinf] = readOptoData(fullfile(DataPath, DataFile), 'channel', channel);
% Fix test info
Dinf = correctTestType(Dinf);
fprintf('Test type: %s\n', Dinf.test.Type);



