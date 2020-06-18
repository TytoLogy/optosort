function spikesBySweep = spikesForAnalysisByUnit(obj, fileNum, ...
																	unitNum, varargin)
%-------------------------------------------------------
% get table of spikes for a specific file, unit and by sweep
% this should be deprecated since it doesn't account for channel
%-------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------

%--------------------------------------
% process options and inputs
%--------------------------------------
% check that file is in range
if ~between(fileNum, 1, obj.Info.nFiles)
	error('requested file %d out of range [1 %d]', ...
							fileNum, obj.Info.nFiles);
elseif ~any(unitNum == obj.listUnits)
	% and that unitnum is valid
	error('unit %d not in Spikes table', unitNum);
end
% check alignment mode
if isempty(varargin)
	% default is no shift
	ALIGN = 'original';
else
	if ~any(strcmpi(varargin{1}, {'original', 'file', 'sweep'}))
		% unknown mode
		error('%s: unknown align option %s', varargin{1})
	else
	ALIGN = lower(varargin{1});
	end
end

%--------------------------------------
% process data
%--------------------------------------
% get valid rows for unit
unit_rows = (obj.Spikes.Unit == unitNum);
% get valid rows for file - this is done by finding the spike
% times that occurred between the specified file's start and end
% times
file_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fileNum)) & ...
						(obj.Spikes.TS <= obj.Info.fileEndTime(fileNum));
% get reduced Spikes table for specified unit and data file
vS = obj.Spikes( (unit_rows & file_rows), :);

% number of sweeps for this file
nsweeps = length(obj.Info.sweepStartTime{fileNum});
% some checks
if nsweeps ~= length(obj.Info.sweepEndTime{fileNum})
	error('mismatch in lengths of sweepStartTime and sweepEndTime');
elseif nsweeps == 0
	error('no sweeps!');	
end
% alignment issues
switch ALIGN
	case 'original'
		% don't modify time stamp values
		alignval = zeros(nsweeps, 1);
	case 'file'
		% adjust timestamp value by first sweep start time
		alignval = obj.Info.sweepStartTime{fileNum}(1) ...
							* ones(nsweeps, 1);
	case 'sweep'
		% adjust timestamp value by each sweep start time
		alignval = obj.Info.sweepStartTime{fileNum};
end

% allocate cell to store spike info for each sweep
spikesBySweep = cell(nsweeps, 1);

% loop through sweeps
for s = 1:nsweeps
	% find the valid time stamps (between sweepStartTimes and
	% sweepEndTimes), store the row indices...
	valid_rows = (vS.TS >= obj.Info.sweepStartTime{fileNum}(s)) ...
					& (vS.TS < obj.Info.sweepEndTime{fileNum}(s));
	% get Spikes that match
	spikesBySweep{s} = vS(valid_rows, :);
	% apply offset correction to timestamps (TS variable)
	spikesBySweep{s}{:, 'TS'} = spikesBySweep{s}{:, 'TS'} - alignval(s);
end

