classdef PLXData
%------------------------------------------------------------------------
% TytoLogy:Experiments:optoproc:PLXDdata class
%------------------------------------------------------------------------
% Class to encapsulate data and code that uses the readPLXFileC()
% function by Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%
%	readPLXFileC:
% 	Author: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% 	Copyright (c) 2012-2013
% 	Last Modified: 2016-10-07 22:29:18
% 	Revision: 4886
%------------------------------------------------------------------------
% Created: 22 April 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: 
%------------------------------------------------------------------------

	properties
		pathname = '';
		filename = '';
		P
	end
	properties (Dependent)
		plxfile
	end
	
	methods

		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		function obj = PLXData(varargin)
		% inputs:
		%	PLXFile		path and filename	
			% if no inputs, we're done
			if isempty(varargin)
				return;
			end
			% get parts of the provided file
			[tmpPath, tmpFile, tmpExt] = fileparts(varargin{1});
			% check if it is a .plx file
			if ~strcmpi(tmpExt, '.plx')
				% ...if not, warn user
				warning('PLXData: file %s does not have .plx extension!', ...
																					varargin{1});
			end
			% assign path and filename.
			obj.pathname = tmpPath;
			obj.filename = [tmpFile tmpExt];
		end
		
		% load data from file, without continuous data, store in P property
		function obj = loadPLX_nocontinuous(obj)
			plxfile = obj.plxfile;
			if isempty(plxfile)
				error('PLXData: filename not defined');
			end
			% read in data using readPLXFileC
			obj.P = readPLXFileC(plxfile, 'all', 'nocontinuous');
		end
		
		%-------------------------------------------------
		%-------------------------------------------------
		% get/set access for dependent properties
		%-------------------------------------------------
		%-------------------------------------------------
		% returns test.Type
		function val = get.plxfile(obj)
			% check to make sure file is specified
			if isempty(obj.filename)
				warning('PLXData: empty filename');
				val = '';
				return
			end
			% create full path and file
			val = fullfile(obj.pathname, obj.filename);
			if ~exist(val, 'file')
				warning('PLXData: file %s not found', val);
				return
			end
		end			
		
	end
	
end