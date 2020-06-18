function spikesBySweep = spikesForAnalysis(obj, fileNum, varargin)
%-------------------------------------------------------
% get table of spikes for a specific file aligned to file start,
% sweep, or as-is (original)
%
%	spikesBySweep = obj.spikesForAnalysis(fileNum)
%
%	Input args:
%		fileNum	file index for data 
%-------------------------------------------------------

%--------------------------------------
% set defaults
%--------------------------------------
% default is no shift
ALIGN = 'original';
% set channel and unit to empty
channelNum = [];
unitNum = [];

%--------------------------------------
% process options and inputs
%--------------------------------------			
% process options
argI = 1;
while argI <= length(varargin)
	switch upper(varargin{argI})
		case 'ALIGN'
			% check alignment mode
			if ~any(strcmpi(varargin{argI+1}, ...
							{'original', 'file', 'sweep'}))
				% unknown mode
				error(['SpikeData.spikesForAnalysis:' ...
							'unknown align option %s'], varargin{argI+1});
			else
				ALIGN = lower(varargin{argI+1});
			end
			argI = argI + 2;
		case {'CHANNEL', 'CHAN', 'C'}
			% user specified channel option, so get desired list
			channelNum = varargin{argI + 1};
			fprintf('SpikeData.spikesForAnalysis: Channel %d\n', ...
							channelNum);
			argI = argI + 2;
		case {'UNIT', 'UN', 'U'}
			% user specified unit(s) so get them from input
			unitNum = varargin{argI + 1};
			fprintf('SpikeData.spikesForAnalysis: Unit %d\n', ...
							unitNum);
			argI = argI + 2;
		otherwise
			% unknown option provided by user
			error(['SpikeData.spikesForAnalysis:' ...
							'unknown option %s'], varargin{argI});
	end
end

%--------------------------------------
% select valid timestamps/table entries
%--------------------------------------
% check that file is in range
if ~between(fileNum, 1, obj.Info.nFiles)
	error('requested file %d out of range [1 %d]', ...
							fileNum, obj.Info.nFiles);
end

%--------------------------------------
% get spikes table for file, channel, unit combination
%--------------------------------------
vS = obj.selectSpikes(fileNum, channelNum, unitNum);

%--------------------------------------
% process sweeps
%--------------------------------------
% number of sweeps for this file
nsweeps = length(obj.Info.sweepStartTime{fileNum});
% some checks
if nsweeps ~= length(obj.Info.sweepEndTime{fileNum})
	error('mismatch in lengths of sweepStartTime and sweepEndTime');
elseif nsweeps == 0
	error('no sweeps!');	
end
% alignment issues (reference timestamps to start of datafile,
% start of sweep, or leave as is  - referenced to start of merged
% file used for spike sorting)
switch ALIGN
	case 'original'
		% don't modify time stamp values
		fprintf('SpikeData:spikesForAnalysis: no timestamp realignment\n');
		alignval = zeros(nsweeps, 1);
	case 'file'
		% adjust timestamp value by first sweep start time for each
		% file in the overall merged file used for sorting
		fprintf('SpikeData:spikesForAnalysis: align timestamp to file start\n');
		alignval = obj.Info.sweepStartTime{fileNum}(1) ...
							* ones(nsweeps, 1);
	case 'sweep'
		% adjust timestamp value by each sweep start time - this is
		% most useful when doing things like analysis and plots of
		% data 
		fprintf('SpikeData:spikesForAnalysis: align timestamp to sweep start\n');
		alignval = obj.Info.sweepStartTime{fileNum};
end

% allocate cell to store spike info for each sweep
spikesBySweep = cell(nsweeps, 1);

% loop through sweeps
for s = 1:nsweeps
	% find the valid time stamps (between sweepStartTimes and
	% sweepEndTimes)
	valid_rows = (vS.TS >= obj.Info.sweepStartTime{fileNum}(s)) ...
					& (vS.TS < obj.Info.sweepEndTime{fileNum}(s));
	spikesBySweep{s} = vS(valid_rows, :);
	% apply offset correction
	spikesBySweep{s}{:, 'TS'} = spikesBySweep{s}{:, 'TS'} - alignval(s);
end
		
