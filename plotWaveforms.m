function H = plotWaveformsFromSpikeTable(varargin)
%------------------------------------------------------------------------
% H = plotWaveformsFromSpikeTable(spikeTable, opts)
%------------------------------------------------------------------------
% TytoLogy:optosort
%--------------------------------------------------------------------------
% THIS MAY NOT  BE NEEDED!!!!!!
%------------------------------------------------------------------------
% Given input table of spike data, plot overlaid waveforms and mean
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 7 July, 2020 (SJS)
%	 
% Revisions:
%
%------------------------------------------------------------------------


if isempty(varargin)
	error('plotWaveFormsFromSpikeTable: need spike table!');
else
	S = varargin{1};
end

% (1) vertically concatenate data stored in cell array of sweeps in
% st.spiketable into a single table
tmpT = vertcat(S.spiketable{:});

% (2) extract just the wave field 
tmpwav = tmpT.Wave;

% (3) plot overlaid waveforms

% plot in new figure
hF = figure;

% need time vector for x-axis
[~, nBins] = size(tmpwav);
t_ms = (1000/S.Info.Fs) * (0:(nBins - 1));

% plot waveforms, with mean and legend
% need tmpwav to be in column format - time in rows, indv. units by column
% so send the function the transpose of the tmpwav matrix.
H = plot_spike_waveforms(t_ms, tmpwav', 'MEAN', true, 'LEGEND', true);

% add title to plot
% create title string with 2 rows:
%	filename (either from S struct or S.Info.FileInfo{findx}.F.file
%	channel and unit
tstr = {	S.fileName, ...
			sprintf('Channel %d Unit %d', channel, unit)};
title(tstr, 'Interpreter', 'none');	
