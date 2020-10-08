function varargout = plotAllData(obj, channel, unit, varargin)
%------------------------------------------------------------------------
% [plot handles] = SpikeData.plotAllDataUnitWaveforms(channel, unit, psthBinWidth)
%------------------------------------------------------------------------
% TytoLogy:opto:optoproc:SpikeData class
%------------------------------------------------------------------------
% 
% Plots: Data, waveforms, etc for selected channel and unit
% 
%------------------------------------------------------------------------
% Input Arguments:
%  channel	channel number of data to plot
%  unit		unit number (for given channel) of data to plot
%
%	Options
%	 psthBinWidth	size of bins for computing PSTH, in milliseconds
%                 default: 5 ms
%   'savePlots'		if provided, plots will be saved (no argument needed)
%	 'saveFormat'	format for saving plots
% 			Save options:
% 				'PNG'		.png (default)
% 				'PDF'		.pdf
% 				'FIG'		.fig (saves MATLAB figure)
%	 'closePlots'  if true, close plots after saving (if not saving plots
%					   will remain open); default is false
%
% Output Arguments:
% 	 H		array of plot handles
% 
%------------------------------------------------------------------------
% Examples
% 	spikedataobj.plotAllData(5, 3)
% 		plots data for channel 5, unit 3. default bin size, 
%		no saving of data
% 	
% 	spikedataobj.plotAllData(4, 2, 10)
% 		plots channel 4, unit 2 with psth bin size of 10 ms
% 
% 	spikedataobj.plotAllData(1, 1, 'savePlots')
% 		plots channel 1, unit 1, and saves plots in current working
% 		directory as default format (.png)
% 		
% 	spikedataobj.plotAllData(6, 4, 'savePlots', 'plotPath', 'E:\Data')
% 		plots channel 6 unit 4 data, saves plots as .png (default) 
% 		in the E:\Data directory
%------------------------------------------------------------------------
% See also: SpikeData class
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 7 July, 2020 (SJS)
%	 
% Revisions:
%	13 Jul 202 (SJS): added docs,  
%  7 Oct 2020 (SJS): added closePlots option
%  8 Oct 2020 (SJS): added plotSnips option, maxSnipsToPlot option
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% define some things
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% default psth bin size is 5 ms
psthBinWidth = 5;
% save plots? default is no
savePlots = false;
% default plot path
plotPath = pwd;
% default plot format
saveFormat = 'PNG';
% default closePlots
closePlots = false;
% default plotSnips
plotSnips = true;
% default max # of snippets to plot
maxSnipsToPlot = 5000;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% parse input args
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% check # of var inputs
nvararg = length(varargin);

% loop through var list of inputs
argN = 1;
while argN <= nvararg

	% if this is the first varargin element and if it is a number, assume
	% that this is psthbinwidth value in milliseconds
	if (argN == 1) && isnumeric(varargin{1})
		% this should be psth Bin size
		psthBinWidth = varargin{1};
		argN = argN + 1;
	else
		% otherwise, check the option and act accordingly
		switch upper(varargin{argN})
			case 'SAVEPLOTS'
				% set savePlots to true
				savePlots = true;
				argN = argN + 1;
			
			case 'PLOTPATH'
				% set plot path to input value
				plotPath = varargin{argN + 1};
				argN = argN + 2;
				
			case 'SAVEFORMAT'
				% if format matches accepted types, assign it to saveFormat
				% variable (as upper case to avoid confusion). 
				if any(strcmpi(varargin{argN + 1}, {'PNG', 'PDF', 'FIG'}))
					saveFormat = upper(varargin{argN + 1});
					argN = argN + 2;
				else
					% if format is invalid, throw error
					error('SpikeData.plotAllData: unknown save format %s', ...
								varargin{argN + 1});
				end
				
			case 'CLOSEPLOTS'
				% set closePlots to true
				closePlots = true;
				argN = argN + 2;

			case 'PSTHBINWIDTH'
				val = varargin{argN + 1};
				if isnumeric(val)
					if val > 0
						% this should be psth Bin size
						psthBinWidth = varargin{1};
					else
						error('%s: psthBinWidth must be greater than 0', ...
									mfilename);
					end
				else
					error('%s: psthBinWidth must be a number', ...
									mfilename);
				end
				argN = argN + 2;
				
			case {'PLOTSNIPS', 'PLOTWAVFORM', 'PLOTSPIKES'}
				val = varargin{argN + 1};
				if all([~isnumeric(val) ~islogical(val)])
					error('%s: arg to plotSnips must be numeric or logical', ...
						      mfilename);
				end
				if val
					plotSnips = true;
				else
					plotSnips = false;
				end
				argN = argN + 2;

			case {'MAXSNIPSTOPLOT', 'MAXSNIPS'}
				val = varargin{argN + 1};
				if ~isnumeric(val)
					error('%s: arg to maxSnipsToPlot must be number', ...
						      mfilename);
				elseif val < 1
					error('%s: arg to maxSnipsToPlot must be >= 1', ...
						      mfilename);
				else
					maxSnipsToPlot = val;
				end
				argN = argN + 2;
				
			otherwise
				error('SpikeData.plotAllData: unknown option %s', ...
							varargin{argN});
		end
	end
end

%------------------------------------------------------------------------
% create output path if necessary
%------------------------------------------------------------------------
if savePlots
	if ~exist(plotPath, 'dir')
		mkdir(plotPath);
	end
	fprintf('%s files will be written to:\n\t%s\n', saveFormat, plotPath);
end

%------------------------------------------------------------------------
% get list of tests
%------------------------------------------------------------------------
testList = cell(obj.Info.nFiles, 1);
for f = 1:obj.Info.nFiles
	testList{f} = obj.Info.FileInfo{f}.testname;
end

%------------------------------------------------------------------------
% loop through tests
%------------------------------------------------------------------------
for findx = 1:obj.Info.nFiles
	sendmsg(sprintf('Plotting data for %s', ...
							obj.Info.FileInfo{findx}.F.file));
						
	%----------------------------------------------
	% get spikes for this file
	%----------------------------------------------
	S = obj.getSpikesByStim(findx, channel, unit);
	% make a local copy of Dinf for this file to make things a little simpler
	Dinf = obj.Info.FileInfo{findx}.Dinf;

	%----------------------------------------------
	% plot waveforms
	%----------------------------------------------
	if plotSnips
		% (1) vertically concatenate data stored in cell array of sweeps in
		% st.spiketable into a single table
		tmpT = vertcat(S.spiketable{:});
		% (2) extract just the wave field 
		tmpwav = tmpT.Wave;
		% (3) plot overlaid waveforms
		% plot in new figure
		% need time vector for x-axis
		[~, nBins] = size(tmpwav);
		t_ms = (1000/obj.Info.Fs) * (0:(nBins - 1));
		% plot waveforms, with mean and legend
		% need tmpwav to be in column format - time in rows, indv. 
		% units by column so send the function the transpose of the 
		% tmpwav matrix.
		hWF = figure;
		hWF_ax = plot_spike_waveforms(t_ms, tmpwav', ...
										'MEAN', true, 'LEGEND', true); %#ok<NASGU>
		% add title to plot
		% create title string with 2 rows:
		%	filename (either from S struct or S.Info.FileInfo{findx}.F.file
		%	channel and unit
		tstr = {	S.fileName, sprintf('Channel %d Unit %d', channel, unit)};
		title(tstr, 'Interpreter', 'none');

		% set figure filename - use the base from the FreqTuningInfo.F object
		set(hWF, 'Name', sprintf('%s_Snips_Ch%d_Un%d', ...
						obj.Info.FileInfo{findx}.F.base, S.channel, S.unit));
		drawnow
		% save waveform plot file
		if savePlots
			fprintf('Saving plot:\n  %s...', ...
								fullfile(plotPath, get(gcf, 'Name')));
			save_plot(hWF, saveFormat, plotPath);
			fprintf('\n...done\n');
			if closePlots
				close(hWF);
			end
		end
	end
	
	%----------------------------------------------
	% plot data according to test
	%----------------------------------------------
	switch upper(obj.Info.FileInfo{findx}.testname)
		case {'BBN', 'LEVEL'}
			% set analysis window to [stimulus onset   stimulus offset]
			analysisWindow = Dinf.audio.Delay + [0 Dinf.audio.Duration];
			% compute rate level function
			RLF = computeRLF(S.spiketimes, S.unique_stim, analysisWindow);
			% plot
			hRLF = plotCurveAndCI(RLF, 'mean');
			% create title string with 2 rows:
			%	filename
			%	channel and unit
			tstr = {	S.fileName, ...
						sprintf('Channel %d Unit %d', channel, unit)};
			% add title to plot
			title(tstr, 'Interpreter', 'none');
			% set filename - use the base from the FreqTuningInfo.F object
			set(hRLF, 'Name', sprintf('%s_Ch%d_Un%d', ...
							obj.Info.FileInfo{findx}.F.base, S.channel, S.unit));
			drawnow
			% save plot
			if savePlots
				fprintf('Saving plot:\n  %s...', ...
									fullfile(plotPath, get(hRLF, 'Name')));
				save_plot(hRLF, saveFormat, plotPath);
				fprintf('\n...done\n');
				if closePlots
					close(hRLF);
				end
			end
			
		case 'FREQ_TUNING'
			% set analysis window to [stimulus onset   stimulus offset]
			analysisWindow = Dinf.audio.Delay + [0 Dinf.audio.Duration];
			% compute rate level function
			FTC = computeFTC(S.spiketimes, S.unique_stim, analysisWindow);
			% plot
			hFTC = plotCurveAndCI(FTC, 'median');
			% create title string with 2 rows:
			%	filename
			%	channel and unit
			tstr = {	S.fileName, ...
						sprintf('Channel %d Unit %d', channel, unit)};
			% add title to plot
			title(tstr, 'Interpreter', 'none');
			% set filename - use the base from the FreqTuningInfo.F object
			set(hFTC, 'Name', ...
							sprintf('%s_Ch%d_Un%d', ...
												obj.Info.FileInfo{findx}.F.base, ...
												S.channel, S.unit));	
			drawnow
			% save plot
			if savePlots
				fprintf('Saving plot:\n  %s...', ...
									fullfile(plotPath, get(hFTC, 'Name')));
				save_plot(hFTC, saveFormat, plotPath);
				fprintf('\n...done\n');
				if closePlots
					close(hFTC);
				end
			end

		
		case 'WAV'
			% plotPSTH is a method in the WAVInfo class, stored at the findx
			% element within the FileInfo cell array. 'LEVEL' tells the method
			% to separate plots according to stimulust level
			hWAV = obj.Info.FileInfo{findx}.plotPSTH(S, psthBinWidth, ...
																	'LEVEL');

			% rename plots with filename, level, channel and unit
			for h = 1:length(hWAV)
				% build on original plot name, since this has the level
				% information already in it
				origname = get(hWAV{h}, 'Name');
				set(hWAV{h}, 'Name', ...
						sprintf('%s_Ch%d_Un%d', origname, S.channel, S.unit));
				drawnow
				% save plot
				if savePlots
					fprintf('Saving plot:\n  %s...', ...
										fullfile(plotPath, get(hWAV{h}, 'Name')));
					save_plot(hWAV{h}, saveFormat, plotPath);
					fprintf('\n...done\n');
					if closePlots
						close(hWAV{h});
					end
				end
			end
			
		case 'FRA'
			% window for spike count
			frawin = Dinf.audio.Delay + [0  Dinf.audio.Duration];
			% calculate FRA stored in struct FRA
			FRA = computeFRA(S.spiketimes, S.unique_stim{1}, ...
														S.unique_stim{2}, frawin);
			% set fname in FRA struct to data file name
			FRA.fname = S.fileName;
			hFRA = plotFRA(FRA, 'dB');	
			% set filename - use the base from the FreqTuningInfo.F object
			set(hFRA, 'Name', sprintf('%s_Ch%d_Un%d', ...
							obj.Info.FileInfo{findx}.F.base, S.channel, S.unit));
			drawnow
			% save plot
			if savePlots
				fprintf('Saving plot:\n  %s...', ...
									fullfile(plotPath, get(hFRA, 'Name')));
				save_plot(hFRA, saveFormat, plotPath);
				fprintf('\n...done\n');
				if closePlots
					close(hFRA);
				end
			end
	
		otherwise
			error('plotAllData: unknown test name %s', ...
									obj.Info.FileInfo{f}.testname)
	end
	

end

%------------------------------------------------------------------------
% eventually need to add all handles to this output
%------------------------------------------------------------------------
if nargout
	if ~closePlots
		varargout{1} = hWF;
	else
		varargout{1} = [];
	end
end