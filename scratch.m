% testing channel, unit sel for SpikeData:spikesForAnalysis

%------------------------------------------------------------------------
% sorted data locations
%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_BBN_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_BBN.nex';
plxFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';
obj = import_from_plexon(fullfile(sortedPath, plxFile), ...
									fullfile(nexPath, nexInfoFile));

%--------------------------------------
%% fake inputs
%--------------------------------------
fileNum =1;
V = {'channel', 4, 'unit', 1}


%--------------------------------------
% set defaults
%--------------------------------------
% default is no shift
ALIGN = 'original';
% set channel and unit to empty
channelNum = [];
unitNum = [];

%--------------------------------------
%% process options and inputs
%--------------------------------------			
% process options
argI = 1;
while argI <= length(V)
	switch upper(V{argI})
		case 'ALIGN'
			% check alignment mode
			if ~any(strcmpi(V{argI+1}, ...
							{'original', 'file', 'sweep'}))
				% unknown mode
				error(['SpikeData.spikesForAnalysis:' ...
							'unknown align option %s'], V{argI+1});
			else
				ALIGN = lower(V{argI+1});
			end
			argI = argI + 2;
		case {'CHANNEL', 'CHAN'}
			% user specified channel option, so get desired list
			channelNum = V{argI + 1};
			fprintf('SpikeData.spikesForAnalysis: Channel %d\n', ...
							channelNum);
			argI = argI + 2;
		case {'UNIT', 'UN'}
			% user specified unit(s) so get them from input
			unitNum = V{argI + 1};
			fprintf('SpikeData.spikesForAnalysis: Unit %d\n', ...
							unitNum);
			argI = argI + 2;
		otherwise
			% unknown option provided by user
			error(['SpikeData.spikesForAnalysis:' ...
							'unknown option %s'], V{argI});
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


% portion of table for this file that contains these channels and units
vS = obj.selectSpikes(fileNum, channelNum, unitNum)

% create temp table of data for desired file
tmpS = obj.spikesForFile(fileNum);

%--------------------------------------
%% check channels and units
%--------------------------------------
% if channelNum is empty, use all channels and units
if isempty(channelNum)
	fprintf(['SpikeData.spikesForAnalysis:' ...
						'using all channels and units\n']);
	channel_rows = true(size(tmpS.Channel));
% 	unit_rows = true(size(tmpS.Channel)); 
	% explicitly set unitNNum to empty
	unitNum = [];
else
	% check channel provided
	cchk = obj.check_channels(channelNum);
	if any(cchk == -1)
		fprintf('SpikeData.spikesForAnalysis: invalid channelNum %d\n', ...
						channelNum(cchk == -1));
		error('SpikeData.spikesForAnalysis: invalid channel');
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
	fprintf(['SpikeData.spikesForAnalysis:' ...
						'using all units for channel %d\n'], channelNum);
	% set unit_rows to ones, size of tmpS.Unit
	unit_rows = true(size(tmpS.Unit));
else
	% if unit num is specified and more than one channel is provided, 
	% throw error if unitNum is not a cell
	if (length(channelNum) > 1) && ~iscell(unitNum)
		error(['SpikeData.spikesForAnalysis: if multiple channels are ' ...
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


%% reduce table to valid channel and unit
vS = tmpS( (channel_rows & unit_rows), :);

%% check units for channel(s)

channels = [4 5];
units = {[0 1], [1 2]};

[uchk, uchklist] = obj.check_units(channels, units)
if any(uchk == false)
	error('Requested unit does not exist on channel');
end
%% code for check_units method for SpikeData (moved into spikeData.m code


if isempty(units)
	error('SpikeData.check_units: list of units to check is empty');
end

if length(channels) > 1
	% if more than one channel is specified, units needs to be a cell
	if ~iscell(units)
		% throw error
		error(['SpikeData.check_units: if multiple channels are ' ...
					'specified, cell vector of units to match must be' ...
					'provided to spikesForAnalysis']);
	else
		% no error - assign units cell to uCell;
		uCell = units;
	end
elseif ~iscell(units)
	% if length channels is 1, put units into a cell to simplify search/check
	% code
	uCell = {units};
end

% check for length mismatch in channels and uCell
if length(channels) ~= length(uCell)
	error(['SpikeData.check_units: mismatch in length of ' ...
					'channels and units']);
end

% allocate output
uList = cell(size(uCell));
val = true(size(uCell));

% get list of units and the corresponding list of channels for the units
[allUnits, allChannels] = obj.listUnits;
% loop through channels
for c = 1:length(channels)
	% get list of units for this channel
	units_for_chan = allUnits{channels(c) == allChannels};
	% check if input units for this channel match
	uchk = ismember(uCell{c}, units_for_chan);
	% add list of zeros to output uList for this channel
	uList{c} = zeros(size(uchk));
	if ~all(uchk)
		% output for val for this channel is false
		val(c) = false;
		for u = 1:length(uCell{c})
			% assign -1 to bad units
			if uchk(u) == false
				uList{c}(u) = -1;
				warning(['SpikeData.check_units: unit %d not found ' ...
						'for channel %d'], uCell{c}(u), channels(c))
			end
		end
	end
end



