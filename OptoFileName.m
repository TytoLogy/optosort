classdef OptoFileName
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%--------------------------------------------------------------------------
	properties
		file
		path
		testdata_file
		base
		animal
		datecode
		unit
		penetration
		depth
		other
	end
	
	methods
		function obj = OptoFileName(varargin)
			if isempty(varargin)
				return
			end
			[p, f, ext] = fileparts(varargin{1});
			obj.path = p;
			obj.file = [f ext];
			f = parseOptoFileName(obj.file);
			obj.base = f.base;
			obj.animal = f.animal;
			obj.datecode = f.datecode;
			obj.unit = f.unit;
			obj.penetration = f.penetration;
			obj.depth = f.depth;
			obj.other = f.other;
		end
	end
end

function F = parseOptoFileName(filename)
	% Input Arguments:
	% 	filename		data file from opto program
	%					e.g., 1372_20191126_03_01_1500_FREQ_TUNING.dat
	% 
	% Output Arguments:
	% 	F	struct with fields:
	% 		animal: '1372'
	% 		datecode: '20191126'
	% 		unit: '03'
	% 		penetration: '01'
	% 		depth: '1500'
	% 		other: 'FREQ_TUNING'
	if isempty(filename)
		F = [];
		return
	end

	%---------------------------------------------------------------------
	% get info from filename
	%---------------------------------------------------------------------
	% only need base filename (no path or ext)
	[~, fname] = fileparts(filename);
	F.base = fname;
	% locate underscores
	usc = find(fname == '_');
	% last underscore index
	endusc = usc - 1;
	% first underscore index
	startusc = usc + 1;
	% animal #
	F.animal = fname(1:endusc(1));
	% date (<year><month><date> => YYYYMMDD, e.g., 20170401 for April 1, 2017)
	F.datecode = fname(startusc(1):endusc(2));
	% unit number
	F.unit = fname(startusc(2):endusc(3));
	% penetration number
	F.penetration = fname(startusc(3):endusc(4));
	% recording depth
	F.depth = fname(startusc(4):endusc(5));
	% this should be test name
	if startusc(4) == startusc(end)
		F.other = fname(startusc(end):end);
	else
		F.other = fname(startusc(5):end);
	end

end
