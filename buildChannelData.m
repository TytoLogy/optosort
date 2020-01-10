function [cD, startI, endI] = buildChannelData(Channels, BPfilt, D, Dinf)
%------------------------------------------------------------------------
% Output = function_template(Input)
%------------------------------------------------------------------------
% <project>:<subproject>:<function_name>
%------------------------------------------------------------------------
% 
% get data for each channel
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Input		input info
% 
% Output Arguments:
% 	Output	output info
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 January, 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% initialize sweep indices to empty arrays
startI = [];
endI = [];

% initialize cD to a store sweeps for each channel
cD = cell(length(Channels), Dinf.nread);

% loop through channels
for c = 1:length(Channels)
	channel = Channels(c);
	% loop through each sweep
	for s = 1:Dinf.nread
		% assign channel sweep data to cD after filtering
		% need to transpose to row vector
		cD{c, s} = D{s}.datatrace(:, channel)';
		% remove initial offset
		cD{c, s} = cD{c, s} - cD{c, s}(1);
		% filter data
		cD{c, s} = filtfilt(BPfilt.b, BPfilt.a, ...
									sin2array(cD{c, s}, ...
										Dinf.indev.Fs, BPfilt.ramp));
		% plot raw and filtered data
		plot(D{s}.datatrace(:, channel)', 'k');
		hold on
			plot(cD{1, s}, 'b');
		hold off
		drawnow
	end

	% build list of sweep start and end indices (in units of samples)
	if isempty(startI)
		% initialize startI and endI to store stop/start
		% locations
		startI = zeros(1, Dinf.nread);
		endI = zeros(1, Dinf.nread);
		% store index points
		if s ~= 1
			startI(s) = endI(s-1) + 1;
		else
			startI(s) = 1;
		end
		endI(s) = startI(s) + length(cD{c, s});
	end
end


