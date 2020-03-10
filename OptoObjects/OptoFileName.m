classdef OptoFileName
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ? 2020 (SJS)
%
% Revisions:
%	3 March 2020 (SJS): altered OptoFileName to deal with plexon sorted 
%								data without "other" element in name
%	10 Mar 2020 (SJS): added newname method to build new filename with
%							additional string and new extension
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

	properties
		file
		path
		testfile
		base
		animal
		datecode
		unit
		penetration
		depth
		other
	end
	
	methods
		%-------------------------------------------------
		% Constructor
		%-------------------------------------------------
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
		
		%-------------------------------------------------
		% return base name without "other" (test name), e.g.
		%	1372_20191126_03_01_1500
		%-------------------------------------------------
		function str = fileWithoutOther(obj)
			str = [	obj.animal '_' ...
								obj.datecode '_' ...
								obj.unit '_' ...
								obj.penetration '_' ...
								obj.depth ];
		end
		
		%-------------------------------------------------
		% utility to create new file using base
		% e.g., [base]_append_str.ext_str
		%-------------------------------------------------
		function str = newname(obj, append_str, ext_str)
			str = [obj.base '_' append_str '.' ext_str];
		end
	end
end

function F = parseOptoFileName(filename)
	%---------------------------------------------------------------------
	% F = parseOptoFileName(filename)
	%---------------------------------------------------------------------
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

	% get info from filename

	% only need base filename (no path or ext)
	[~, fname] = fileparts(filename);
	F.base = fname;
	% locate underscores
	usc = find(fname == '_');
	if isempty(usc)
		st = dbstack;
		error('%s: Improper Opto File Name form %s', mfilename, st.name);
	end
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
	% need to check if we're at final underscor
	if length(usc) == 4
		% only thing left is rec depth and then . and extension
		F.depth = fname(startusc(4):end);
		% no other
		F.other = '';
	else
		% recording depth
		F.depth = fname(startusc(4):endusc(5));
		% this should be test name
		if startusc(4) == startusc(end)
			F.other = fname(startusc(end):end);
		else
			F.other = fname(startusc(5):end);
		end
	end
end
