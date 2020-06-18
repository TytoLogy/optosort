%-------------------------------------------------------
%-------------------------------------------------------
function H = plotUnitWaveforms(obj, channel, varargin)
%-------------------------------------------------------
% [plot handles] = obj.plotUnitWaveforms([channels], 
% Plot sorted waveforms for each identified unit for a given channel
% This will work for individual channels and either all units for the
% channel (if unit list is not provided) or a specified unit(s)
%-------------------------------------------------------
if length(channel) > 1
	error('SpikeData.plotUnitWaveforms: single channel only');
end
% check channel provided
cchk = obj.check_channels(channel);
if cchk == -1
	error('SpikeData.plotUnitWaveforms: invalid channel %d', ...
						channel);
end
if isempty(varargin)
	% use all units for this channel
	tmp = obj.listUnits(channel);
	unitList = tmp{1};
	nU = length(unitList);
else
	% use list provided by user
	unitList = varargin{1};
	% check units
	nU = length(unitList);
	if nU == 0
		error(['SpikeData.plotUnitWaveforms: ' ...
					'no units for channel %d'], channel);
	else
		% check units
		unitsForChannel = obj.listUnits(channel);
		for u = 1:nU
			if ~any(unitList(u) == unitsForChannel{1})
				error('unit %d not found', unitList(u));
			end
		end
	end
end			
% allocate gobjects array to hold figure handles
H = gobjects(nU, 1);
% loop through units
for u = 1:nU
	fprintf('Plotting unit %d waveforms\n', unitList(u));
	% create figure and store handle in H array
	H(u) = figure;
	% get spike waveforms for this unit
	W = obj.Spikes{obj.Spikes.Unit == unitList(u), 'Wave'};
	if ~isempty(W)
		[~, nBins] = size(W);
		ms = (1000/obj.Info.Fs) * (0:(nBins - 1));
		plot(ms, W', 'k');
	end
	[~, fstr] = fileparts(obj.Info.FileName);
	tstr = { sprintf('File: %s', [fstr '.mat']), ...
				sprintf('Channel %d ', channel), ...
				sprintf('Unit %d ', unitList(u)) };
	title(tstr, 'Interpreter', 'none');
	xlabel('Time (ms)');
% 				grid on
	H(u).Name = sprintf('%s_unit%d', fstr, unitList(u));
end
% make box tight around data
axis tight
