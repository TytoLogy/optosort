test = Dinf.test;
test.audio = Dinf.audio;
test.opto = Dinf.opto;
tdt.outdev = Dinf.outdev;
test.audio.signal.Type = char(Dinf.audio.signal.Type);


[stimcache, stimseq] = opto_buildStimCache(test, tdt, Dinf.caldata);