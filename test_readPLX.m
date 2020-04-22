% script to test/play with reading .plx (plexon) files using readPLXFile
% created by Benjamin Kraus

% readPLXFile help:
% readPLXFileC - A MEX function to read a PLX file (Plexon, Inc.).
% 
% USAGE:
% plx = readPLXFileC(filename, varargin)
% plx = readPLXFileC('help')
% plx = readPLXFileC('version')
% 
% INPUT:
% filename - Name of the PLX file to read.
% varargin - One (or more) of the arguments listed below. Arguments are
% 			  parsed in order, with later arguments overriding earlier
% 			  arguments.
% 
% ARGUMENTS:
% 'help'           - Display this help information
% 'version'        - Display MEX file version information
% 						 If 'version' occurs as the first input argument,
% 						 the revision number is returned as the first (and only) output,
% 						 and the version information is only printed to screen
% 						 if no ouptut is requested.
% 						 If 'version' occurs after the first input argument,
% 						 version information is printed to the screen, but
% 						 otherwise the function behaves as though 'version' was not present.
% 'headers'        - Retrieve only headers (default)
% 						 (implies 'nospikes','noevents','nocontinuous')
% '[no]fullread'   - Scan the entire file (default = 'nofullread')
% 						 ('fullread' is implied if anything other than headers are requested)
% '[no]spikes'     - Retrieve (or not) spike timestamps (default = 'nospikes')
% 						 'nospikes' implies 'nowaves'
% '[no]waves'      - Retrieve (or not) spike waveforms (default = 'nowaves')
% 						 'waves' implies 'spikes'
% '[not]units'     - Must be followed by a list of units to (not) retrieve
% 						 0 = unsorted, 1 = unit 'a', 2 = unit 'b', etc.
% '[no]events'     - Retrieve (or not) event data (default = 'noevents')
% '[no]continuous' - Retrieve (or not) continuous data (default = 'no')
% 'all'            - Read the entire file
% 						 (implies 'spikes','waves','events','continuous')
% 'range'          - Time range of data to retrieve
% 'start'          - Start of time range of data to retrieve
% 'stop'           - End of time range of data to retrieve
% 'first'          - First data sample to retrieve
% 'num'            - Number of data samples to retieve
% 'last'           - Last data sample to retrieve
% 
% SELECTING CHANNELS:
% 'spikes','waves','events', and/or 'continuous' can be followed by a
% numerical array, which is then parsed to determine which channels to
% retrieve. An empty array implies 'no'. If the array is missing,
% then all channels are retrieved.
% 
% OUTPUT:
% plx - A structure containing the PLX file data.

%% path to readPXFileC

if ~exist('readPLXFileC')
	fprintf('adding readPLXFile to path\n');
	addpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/Plexon/readPLXFileC');
end

%% define input files

% path to plx file(s)
PLXPath = '/Volumes/Lexar/Work/DATA/Plexon-ICdata';
% plx file to open
PLXFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';
fname = fullfile(PLXPath, PLXFile);

%% read in everything
p = readPLXFileC(fname, 'fullread', 'all');

%%
p2 = readPLXFileC(fname, 'all', 'nocontinuous')

%% test object

pObj = PLXData(fname);

%%
