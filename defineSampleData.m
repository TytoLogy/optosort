%------------------------------------------------------------------------
%------------------------------------------------------------------------
function varargout = defineSampleData(varargin)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% defineSampleData.m
%------------------------------------------------------------------------
% defines path to sample data for testing
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 18 May, 2020 (extracted from export_for_plexon) (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


Channels = [];

% check inputs 
if nargin == 0
	% do something to get data files from user !!!!
	DataPath = {};
	DataFile = {};
	TestFile = {};
	Channels = [];
elseif nargin == 3
	% assign values to DataPath, DataFile and TestFile
	DataPath = varargin{1};
	DataFile = varargin{2};
	TestFile = varargin{3};
else
	error('%s->defineSampleData: invalid inputs', mfilename)
end

% loop through # of data files, create file objects
for f = 1:length(DataFile)
	if length(DataPath) == 1
		% only 1 element in DataPath so assume all data files are on this
		% path
		F(f) = OptoFileName(fullfile(DataPath{1}, DataFile{f})); %#ok<AGROW>
	else
		F(f) = OptoFileName(fullfile(DataPath{f}, DataFile{f})); %#ok<AGROW>
	end
	F(f).testfile = TestFile{f};  %#ok<AGROW>
end	

% assign outputs
varargout{1} = F;
if nargout == 2
	varargout{2} = Channels;
end
