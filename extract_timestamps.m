function out = extract_timestamps(spikeTableCell)
%------------------------------------------------------------------------
% out = extract_timestamps(spikeTableCell)
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%------------------------------------------------------------------------
% Given a spike table, convert to a cell array of spike times for
% processing/plotting
%------------------------------------------------------------------------
% spikeTableCell		sorted spikes in table
% 				Table Variable Names:
% 					Channel		AD channel
% 					Unit			Unit ID (for given channel! note that units might
% 									 not have unique IDs across channels)
% 					TS				timestamp (seconds)
% 					PCA			PCA values (not valid for data imported 
% 									 directly from plx file
% 					Wave			Wave snippet data					
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 6 May 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
%------------------------------------------------------------------------
nT = length(spikeTableCell);
out = cell(nT, 1);
for n = 1:nT
	out{n} = spikeTableCell{n}.TS;
end
	
