function varargout = plot_snippets(SD, SpikeWindow, Fs)
%------------------------------------------------------------------------
% plot_snippets(SD)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optoproc
%--------------------------------------------------------------------------
% plots snippets (spike waveforms) from struct returned by 
% threshold_opto_data()
%
%------------------------------------------------------------------------
% Inputs:
%	SD		struct with fields:
% 						ts  cell array of spiketimes (in milliseconds)
% 						snips		cell array of spike waveform snippets
%	SpikeWindow	[preTS postTS] time, in milliseconds
%	Fs		Sampling rate for spike data
% Outputs:
%------------------------------------------------------------------------
% See Also: opto, threshold_opto_data
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 March 2020 (SJS)
%
% Revisions:
%--------------------------------------------------------------------------

% get # of stimuli
nstim = length(SD.snips);

if nstim == 0
	error('%s: no stimuli (length of snips == 0)', mfilename);
end

% Plot the snippets

% loop through stimuli (elements of SD)
for n = 1:nstim
	% get the # of reps for this stimulus
	nreps = length(SD.snips{n});
	% and loop through them
	for r = 1:nreps
		% if first plot, create new plot
		if n == 1
			% need to turn the snippets into a column-oriented array using
			% transpose
			plot(SD.snips{n}{r}', 'k');
		else
			% hold and plot
			hold on
				plot(SD.snips{n}{r}', 'k');
			hold off
		end
	end
end

% turn on grid
grid on

% offset x axis to account for spike window pre-detect time
offset_bins = ms2bin(SpikeWindow(1), Fs);

% get current tick marks (in samples)
xtbins = get(gca, 'XTick');
% apply offset
xtms = bin2ms( xtbins - offset_bins, Fs);
% re-do tick labels in milliseconds
xtl = get(gca, 'XTickLabel');
for n = 1:length(xtl)
	xtl{n} = sprintf('%.1f', xtms(n));
end
set(gca, 'XTickLabel', xtl);
xlabel('ms');

varargout{1} = gca;
