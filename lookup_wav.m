function wavdata = lookup_wav(wavname, wavs)
% wavs is struct array of wav information from opto program
	
	% create cell list of wav names
	wavnames = {wavs.Filename};
	
	% use name to locate info
	wavdata = wavs(strcmpi(wavname, wavnames));	