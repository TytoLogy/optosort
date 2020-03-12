classdef TestSuperClass
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 12 March 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: how to handle multiple channels: answer: Traces will hold them
%------------------------------------------------------------------------

	properties
		% Info				Info object (CurveInfo or subclass)
		% Traces				cell array of spike data
		% Spikes				spikes, snippets, threshold, channel info, cell array
		% 
		Info
		Traces
		Spikes
	end
% 	properties (Dependent)
% 		Channels
% 		nChannels
% 	end
	methods
		function obj = TestSuperClass(varargin)
			% constructor needs 0 inputs or Dinf, Traces, Spikes
			if length(varargin) ~= 2
				return;
			end
			% create CurveInfo object, store in Info
			obj.Info = CurveInfo(varargin{1});
			obj.Traces = varargin{2};
			% see if spikes were provided
			if length(varargin) == 3
				Spikes = varargin{3}; %#ok<NASGU>
			end
		end
		
		%-------------------------------------------------------
		% return a list of channels in the Spikes array
		%-------------------------------------------------------
		function [Channels, nChan] = listChannels(obj)
			if isempty(obj.Spikes)
				warning('No Spikes loaded')
				Channels = [];
				nChan = 0;
			else
				% need to create array of Spikes structs
				tmp = [obj.Spikes{:}];
				% then get array of channels
				Channels = [tmp.Channel];
				% and count them
				nChan = length(Channels);
			end
		end
		
		%-------------------------------------------------------
		% Plot sorted waveforms for each identified unit
		%-------------------------------------------------------
		function varargout = plotWaveformsForChannel(obj, varargin)
			if ~ obj.hasSpikes
				warning('No Spikes loaded!')
				if nargout
					varargout{1} = [];
				end
				return
			end
 			if isempty(varargin)
				% plot all channels
				channel_to_plot = obj.listChannels;
			else
				channel_to_plot = varargin{1};
			end
			% look for channels
			% first, get channels present in Spikes
			[ChannelsLoaded, ~] = obj.listChannels;
			% then look for indices to the desired channels
			H = zeros(length(channel_to_plot), 1);
			for c = 1:length(channel_to_plot)
				chanID = find(channel_to_plot(c) == ChannelsLoaded);
				if chanID == 0 
					error('channel %d not found in Spikes', channel_to_plot(c));
				end
				figure
				H(c) = plot_snippets(obj.Spikes{chanID}, ...
						obj.Spikes{chanID}.tset.SpikeWindow, obj.Info.ADFs);
				title(	{obj.Info.F.file, ...
								sprintf('Channel %d', channel_to_plot(c))}, ...
							'Interpreter', 'none');
				set(gcf, 'Name', [obj.Info.F.base '-snips-Ch' ...
									num2str(channel_to_plot(c))])
			end
			if nargout
				varargout{1} = H;
			end
		end
		

		function val = hasSpikes(obj)
		%-------------------------------------------------------
		% returns true if Spikes struct has been assigned, false if not
		%-------------------------------------------------------
			if isempty(obj.Spikes)
				val = false;
			else
				val = true;
			end
		end
		
		
		
	end
end
