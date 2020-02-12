classdef SpikeInfo
	properties
		Fs
		sweepStartBin = {};
		sweepEndBin = {};
		fileStartBin
		fileEndBin
	end	
	properties (Dependent)
		fileStartTime
		fileEndTime
		sweepStartTime
		sweepEndTime
	end
	
end
	
		