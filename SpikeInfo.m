classdef SpikeInfo
	properties
		% FileName		nex file name (include path)
		% InfoFileName	nexinfo filename (incluedes path)
		% FileData		struct array of information for each file
		%array of OptoFileName objects holding .dat file info
		% Fs				samples/sec for data
		% sweepStartBin	{nfiles, 1} holding arrays of sample for start of
		%						each data recording sweep
		% sweepEndBin		{nfiles, 1} holding arrays of sample for end of
		%						each data recording sweep
		% fileStartBin		[nfiles, 1] sample for start of
		%						each data file within .nex file
		% fileEndBin		[nfiles, 1] sample for end of
		%						each data file within .nex file
		% startBinVector	[total # sweeps, 1] sample for start of each sweep
		% endBinVector	[total # sweeps, 1] sample for end of each sweep
		% ADchannel		list of AD channels exported within .nex file
		FileName
		InfoFileName
		FileData
		Fs
		sweepStartBin = {};
		sweepEndBin = {};
		fileStartBin
		fileEndBin
		startBinVector
		endBinVector
		ADchannel
		dataFilter
	end	% END properties (main)
	properties (Dependent)
		fileStartTime
		fileEndTime
		sweepStartTime
		sweepEndTime
		startTimeVector
		endTimeVector
		nFiles
	end	% END properties(Dependent)
	
	methods
		%-------------------------------------------------
		%-------------------------------------------------
		% Constructor
		%-------------------------------------------------
		%-------------------------------------------------
		function obj = SpikeInfo(varargin)
			if isempty(varargin)
				return
			end
			if length(varargin) == 2
				if strcmpi(varargin{1}, 'struct')
					% this might be deprecated/removed in future
					obj = obj.initFromNexInfoStruct(varargin{2});
				elseif strcmpi(varargin{1}, 'file')
					obj = obj.initFromNexInfoFile(varargin{2});
				else
					error('Unknown input type %s', varargin{1});
				end
			else
				error('need both input mode and input value');
			end
		end

		%-------------------------------------------------
		%-------------------------------------------------
		% utility methods
		%-------------------------------------------------
		%-------------------------------------------------

		%-------------------------------------------------
		% initialize from NexInfo struct... remove????
		%-------------------------------------------------
		function obj = initFromNexInfoStruct(obj, nexInfo)
			obj.FileName = nexInfo.NexFileName;
			obj.Fs = nexInfo.Fs;
			obj.sweepStartBin = nexInfo.sweepStartBin;
			obj.sweepEndBin = nexInfo.sweepEndBin;
			obj.fileStartBin = nexInfo.fileStartBin;
			obj.fileEndBin = nexInfo.fileEndBin;
		end
		
		function obj = initFromNexInfoFile(obj, nexInfoFileName)
% 			obj = load(nexInfoFileName, 'nexInfo');
%			obj = load(nexInfoFileName);
			if ~exist(nexInfoFileName, 'file')
				error('nexinfo file %s not found', nexInfoFileName);
			else
				tmpStruct = load(nexInfoFileName, 'nexInfo');
				obj = tmpStruct.nexInfo;
				clear tmpStruct
			end
		end
		
		function val = nChannels(obj)
			val = length(obj.ADchannel);
		end
		
		%-------------------------------------------------
		%-------------------------------------------------
		% access for dependent properties
		%-------------------------------------------------
		%-------------------------------------------------
		
		%-------------------------------------------------
		% # of files merged into nex file
		%-------------------------------------------------
		function val = get.nFiles(obj)
			val = length(obj.FileData);
		end
		%-------------------------------------------------
		% fileStart, end Time computed from bins
		%-------------------------------------------------
		function val = get.fileStartTime(obj)
			val = (1/obj.Fs) * (obj.fileStartBin - 1);
		end
		function val = get.fileEndTime(obj)
			val = (1/obj.Fs) * (obj.fileEndBin - 1);
		end
		
		% convert sweep bin cells to vectors... 
		function val = get.startBinVector(obj)
			val = [obj.sweepStartBin{:}];
		end
		function val = get.endBinVector(obj)
			val = [obj.sweepEndBin{:}];
		end
		%-------------------------------------------------
		% sweep Start, end Time computed from bins
		%-------------------------------------------------
		function val = get.sweepStartTime(obj)
			val = cell(obj.nFiles, 1);
			for f = 1:obj.nFiles
				val{f} = (obj.sweepStartBin{f} - 1) * (1/obj.Fs);
			end
		end
		function val = get.sweepEndTime(obj)
			val = cell(obj.nFiles, 1);
			for f = 1:obj.nFiles
				val{f} = (obj.sweepEndBin{f} - 1) * (1/obj.Fs);
			end
		end
		function val = get.startTimeVector(obj)
			val = (obj.startBinVector - 1) * (1/obj.Fs);
		end
		function val = get.endTimeVector(obj)
			val = (obj.endBinVector - 1) * (1/obj.Fs);
		end
	end	% END methods
	
end	% END classdef
	
		