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
%  9 Jul 2020 (SJS): adding multifile gui select
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


Channels = [];

% check inputs 
if nargin == 0
	% initialize path and file arrays
	DataPath = {};
	DataFile = {};
	
	% get .dat file(s)
	% open uigetfile with multiselect enabled
	% https://www.mathworks.com/matlabcentral/answers/84455-open-many-files-using-uiopen-or-uigetfile
	[tmpfile, tmppath] = uigetfile('*.dat', 'Select opto .dat file(s)', ...
											'MultiSelect', 'on');
	% check if user canceled (DataFile will be 0)
	if isnumeric(tmpfile)
		if tmpfile == 0
			varargout{1} = [];
			return;
		end
	else
		% otherwise assign values
		if ischar(tmpfile)
			DataFile{1} = tmpfile;
		else
			DataFile = tmpfile;
		end
		DataPath{1} = tmppath;
	end
	
	% create test file names and see if they can be found
	nf = length(DataFile);
	TestFile = cell(size(DataFile));
	for f = 1:nf
		% create object
		fobj(f) = OptoFileName(fullfile(DataPath{1}, DataFile{f})); %#ok<AGROW>
		% create test filename
		testname = [fobj(f).base '_testdata.mat'];
		if ~exist(fullfile(DataPath{1}, testname), 'file')
			[TestFile{f}, ~] = uigetfile(testname, 'Select _testdata.mat file');
			if isnumeric(TestFile{f})
				if TestFile{f} == 0
					varargout{1} = [];
					return
				end
			end
		else
			TestFile{f} = testname;
		end
	end
	
	% get list of channels
	s = inputdlg('Channels (1-16) (e.g. [1 2:4 6])');
	if isempty(s)
		varargout{1} = [];
		return
	else
		Channels = str2num(s{:}); %#ok<ST2NM>
	end
	
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
