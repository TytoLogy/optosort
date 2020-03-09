function [varlist, nvars] = get_opto_varlist(cInfo)
%---------------------------------------------------------------------
% Some test-specific things...
%---------------------------------------------------------------------
switch upper(cInfo.testtype)
	case {'FREQ', 'LEVEL'}
		% list of frequencies, and # of freqs tested
		% list of levels, and # of levels tested
		varlist = {cInfo.varied_values};
		nvars = length(varlist);

	case 'FREQ+LEVEL'
		% list of freq, levels
		varlist = cell(2, 1);
		% # of freqs in nvars(1), # of levels in nvars(2)
		nvars = zeros(2, 1);
		tmprange = cInfo.varied_values;
		for v = 1:2
			varlist{v} = unique(tmprange(v, :), 'sorted');
			nvars(v) = length(varlist{v});
		end

	case 'OPTO'
		% not yet implemented
		
	case 'WAVFILE'
		% get list of stimuli (wav file names)
		varlist = cInfo.Dinf.test.wavlist;
		nvars = length(varlist);

	otherwise
		error('%s: unsupported test type %s', mfilename, cInfo.testtype);
end