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
S = import_from_plexon(fullfile(sortedPath, plxFile), ...
									fullfile(nexPath, nexInfoFile));

%%

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
% process options and inputs
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

% create temp table of data for desired file
tmpS = obj.spikesForFile(fileNum);

%--------------------------------------
% check channels and units
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
	% set unit_rows to ones, size of tmpS.Channel
	unit_rows = true(size(tmpS.Channel));
elseif unitNum ~= -1
	% get indices for unit
	unit_rows = tmpS.Unit == unitNum;				
end
% reduce table to valid channel and unit
vS = tmpS( (channel_rows & unit_rows), :);
% clear tmpS
clear tmpS