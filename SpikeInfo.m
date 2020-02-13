classdef SpikeInfo
	properties
		% FileName		nex file name
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
				if strcmpi(varargin{2}, 'nexInfo')
					obj = obj.initFromNexInfo(varargin{1});
				else
					error('Unknown input type %s', varargin{1});
				end
			else
				error('need both input struct and ID string');
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
		function obj = initFromNexInfo(obj, nexInfo)
			obj.FileName = nexInfo.NexFileName;
			obj.Fs = nexInfo.Fs;
			obj.sweepStartBin = nexInfo.sweepStartBin;
			obj.sweepEndBin = nexInfo.sweepEndBin;
			obj.fileStartBin = nexInfo.fileStartBin;
			obj.fileEndBin = nexInfo.fileEndBin;
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
			val = length(obj.DataFiles);
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
	
		