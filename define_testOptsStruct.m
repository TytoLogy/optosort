function testOpts = define_testOptsStruct(Channels, varargin)
% testOpts will be used by buildTiminTestData, test_data_for_export, and
% getTestData to build test data "recordings" for each data sweep

testOpts = repmat( struct( 'nUnits', [], ...
									'FiringRate', [], ...
									'BGRate', [], ...
									'UnitAmp', [], ...
									'BurstOnsetTime', [], ...
									'BurstOffsetTime', [], ...
									'Jitter', [] ...
									), ...
							length(Channels), 1);
if ~isempty(varargin)
	Dinf = varargin{1};
else
	Dinf.audio.Delay = 100;
	Dinf.audio.Duration = 150;
end

for c = 1:length(Channels)
	% set # of units per channel
	testOpts(c).nUnits = c;
	% set the peak firing rate per Channel to be equal to channel value
	testOpts(c).FiringRate = 5*Channels(c);
	% background rate for each unit? (unused for now)
	testOpts(c).BGRate = 0;
	% unit amplitude (peak)
	testOpts(c).UnitAmp = ones(testOpts(c).nUnits, 1);
	for u = 1:testOpts(c).nUnits
		testOpts(c).UnitAmp(u) = u * (1/testOpts(c).nUnits);
	end
	% set one channel to start at onset time
	if c == length(Channels)
		testOpts(c).BurstOnsetTime = Dinf.audio.Delay;
		testOpts(c).BurstOffsetTime = Dinf.audio.Delay + Dinf.audio.Duration;
	else
		testOpts(c).BurstOnsetTime = 0;
		testOpts(c).BurstOffsetTime = 400;
	end
end