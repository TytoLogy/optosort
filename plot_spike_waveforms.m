function varargout = plot_spike_waveforms(tvec, wavematrix, varargin)
%------------------------------------------------------------------------
% [plot handles] =  plot_spike_waveforms(tvec, wavematrix, 
%------------------------------------------------------------------------
% TytoLogy:optosort
%------------------------------------------------------------------------
% 
% Plots: Data, waveforms, etc for selected channel and unit
% 
%------------------------------------------------------------------------
% Input Arguments:
%  tvec			[1, nsamples] time vector (in milliseconds... or ... ?)
%  wavematrix	[nsamples, nwaves] matrix of spike waveforms
%  Optional:
%	 'mean'		plot mean waveform (default: true)
%
% Output Arguments:
% 	 [plot handles]		array of plot handles
% 
%------------------------------------------------------------------------

PLOTMEAN = true;
PLOTLEGEND = true;

nvararg = length(varargin);
argI = 1;
while argI <= nvararg
	switch(upper(varargin{argI}))
		case {'MEAN', 'AVG', 'PLOTMEAN', 'PLOTAVG'}
			if varargin{argI+1} == 1
				PLOTMEAN = true;
			else
				PLOTMEAN = false;
			end
			argI = argI + 2;
		case {'LEGEND', 'PLOTLEGEND'}
			if varargin{argI+1} == 1
				PLOTLEGEND = true;
			else
				PLOTLEGEND = false;
			end
			argI = argI + 2;
			
		otherwise
			error('%s: unknown option %s', mfilename, varargin{argI})
	end
end


% plot all spike waveforms, overlaid in grey, saving handle
aH = plot(tvec, wavematrix', 'Color', 0.5*[1 1 1]);

if PLOTMEAN
	% ... and plot average waveform
	hold on
		mH = plot(tvec, mean(wavematrix, 2), 'b');
	hold off
end

% make axis close to range of data
axis tight

% x axis label
xlabel('time (ms)')

if PLOTLEGEND
	% legend
	legend([aH(1), mH], {'indv.', 'mean'}, ...
									'Box', 'off', 'Location', 'best');
end

if nargout
	if PLOTMEAN
		varargout{1} = {aH, mH};
	else
		varargout{1} = {aH};
	end
end
