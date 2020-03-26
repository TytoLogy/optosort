classdef WAVtestdata < TestSuperClass
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 11 March 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: how to handle multiple channels?????
%------------------------------------------------------------------------

	properties
		% Info				WAVInfo object
		% Traces				cell array of spike data
		% Spikes				spikes, snippets, threshold, channel info, cell array
		% 
% 		Info
% 		Traces
% 		Spikes
	end
% 	properties (Dependent)
% 		Channels
% 		nChannels
% 	end
	methods
		function obj = WAVtestdata(varargin)
			% constructor needs 0 inputs or Dinf, Traces, Spikes
			if length(varargin) ~= 2
				return;
			end
			% check Dinf.test.Name to make sure we were passed WAV / WAVFILE
			% test data
			tmpDinf = varargin{1};
			if ~strcmpi('WAV', tmpDinf.test.Name)
				error('%s: Data (Dinf) are not named WAV (Name = %s)', ...
										mfilename, tmpDinf.test.Name)
			end
			% create WAVInfo object, store in Info
			obj.Info = WAVInfo(varargin{1});
			obj.Traces = varargin{2};
			% see if spikes were provided
			if length(varargin) == 3
				Spikes = varargin{3}; %#ok<NASGU>
			end
		end
		
		%{
		%-------------------------------------------------------
		% return a list of channels in the Spikes array
		%-------------------------------------------------------
		function [Channels, nChan] = listChannels(obj)
			if isempty(obj.Spikes)
				Channels = [];
				nChan = 0;
			else
				tmp = [obj.Spikes{:}];
				Channels = [tmp.Channel];
				nChan = length(Channels);
			end
		end
		
		
		%-------------------------------------------------------
		% Plot sorted waveforms for each identified unit
		%-------------------------------------------------------
		function H = plotWaveformsForChannel(obj, aChannel)
			% note Fs is not defined yet!!!!!
			[Channels, ~] = obj.listChannels;
			chanID = find(Channels == aChannel);
			if isempty(chanID)
				error('%s: channel %d not found in Spikes', mfilename, aChannel);
			end
			H = plot_snippets(obj.Spikes{chanID}, ...
						obj.Spikes{chanID}.tset.SpikeWindow, Fs);
		end
		
		%}
		
	end
end
