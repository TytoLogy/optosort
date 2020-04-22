% script to test indexing into wav test

% get list of stimuli (wav file names)
nwavs = length(Dinf.stimList);
wavlist = cell(nwavs, 1);
stimindex = cell(nwavs, 1);
% loop through 
for w = 1:nwavs
	stype = Dinf.stimList(w).audio.signal.Type;
	if strcmpi(stype, 'null')
		wavlist{w} = 'null';
	elseif strcmpi(stype, 'noise')
		wavlist{w} = 'BBN';
	elseif strcmpi(stype, 'wav')
		[~, wavlist{w}] = ...
			fileparts(Dinf.stimList(w).audio.signal.WavFile);
	else
		error('%s: unknown type %s', mfilename, stype);
	end
	stimindex{w} = find(Dinf.test.stimIndices == w);
end