function [tbl, varargout] = selectSpikes(obj, fileNum, channelNum, unitNum)
%-------------------------------------------------------
%-------------------------------------------------------
%-------------------------------------------------------

% create temp table of data for desired file
tmpS = obj.spikesForFile(fileNum);

% check channels and units
% if channelNum is empty, use all channels and units
if isempty(channelNum)
	fprintf(['SpikeData.selectSpikes:' ...
						'using all channels and units\n']);
	channel_rows = true(size(tmpS.Channel));
% 	unit_rows = true(size(tmpS.Channel)); 
	% explicitly set unitNNum to empty
	unitNum = [];
else
	% check channel provided
	cchk = obj.check_channels(channelNum);
	if any(cchk == -1)
		fprintf('SpikeData.selectSpikes: invalid channelNum %d\n', ...
						channelNum(cchk == -1));
		error('SpikeData.selectSpikes: invalid channel');
	end
	% get indices for channel(s)
	channel_rows = false(size(tmpS.Channel));
	% loop through channels and OR channel_rows with channels
	for c = 1:length(channelNum)
		channel_rows = channel_rows | (tmpS.Channel == channelNum(c));
	end
end

% if unitNum is empty, find all units
if isempty(unitNum)
	fprintf(['SpikeData.selectSpikes:' ...
						'using all units for channel %d\n'], channelNum);
	% set unit_rows to ones, size of tmpS.Unit
	unit_rows = true(size(tmpS.Unit));
else
	% if unit num is specified and more than one channel is provided, 
	% throw error if unitNum is not a cell
	if (length(channelNum) > 1) && ~iscell(unitNum)
		error(['SpikeData.selectSpikes: if multiple channels are ' ...
					'specified, cell vector of units to match must be' ...
					'provided to spikesForAnalysis']);
	else
		% check units
		[uchk, uchklist] = obj.check_units(channelNum, unitNum); %#ok<ASGLU>
		if any(uchk == false)
			error('Requested unit does not exist on channel');
		end
	end

	% build unit_rows
	% initialize unit_rows
	unit_rows = false(size(tmpS.Unit));
	% loop through units, OR unit_rows with each unit_num
	for u = 1:length(unitNum)
		% get indices for unit
		unit_rows = unit_rows | (tmpS.Unit == unitNum);
	end
end

% reduce table to valid channel and unit
tbl = tmpS( (channel_rows & unit_rows), :);

% outputs
varargout{1} = channel_rows;
varargout{2} = unit_rows;

