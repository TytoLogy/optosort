classdef SpikeInfo
	properties
		FileName
		DataFiles
		Fs
		sweepStartBin = {};
		sweepEndBin = {};
		fileStartBin
		fileEndBin
		ADchannels
	end	% END properties
	properties (Dependent)
		fileStartTime
		fileEndTime
		sweepStartTime
		sweepEndTime
		nFiles
	end	% END properties(Dependent)
	
	methods
		
		function obj = SpikeInfo(varargin)
			if isempty(varargin)
				return
			end
			if length(varargin) == 2
				if strcmpi(varargin{2}, 'nexInfo')
					obj = obj.initFromNexInfo(varargin{1});
				else
					error('Unknown input type %s', varargin{1});
				end
			else
				error('need both input struct and ID string');
			end
		end
		
		function obj = initFromNexInfo(obj, nexInfo)
			obj.FileName = nexInfo.NexFileName;
			obj.Fs = nexInfo.Fs;
			obj.sweepStartBin = nexInfo.sweepStartBin;
			obj.sweepEndBin = nexInfo.sweepEndBin;
			obj.fileStartBin = nexInfo.fileStartBin;
			obj.fileEndBin = nexInfo.fileEndBin;
		end
		
		function val = get.fileStartTime(obj)
			val = (1/obj.Fs) * (obj.fileStartBin - 1);
		end
		function val = get.fileEndTime(obj)
			val = (1/obj.Fs) * (obj.fileEndBin - 1);
		end
	end	% END methods
	
end	% END classdef
	
		