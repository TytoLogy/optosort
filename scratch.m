[plxfile, plxpath] = uigetfile({	'*.plx', 'Plexon file (*.plx)'; ...
											'*.*', 'All Files (*.*)' }, ...
											'Select .plx file');
% return if cancelled
if plxfile == 0
	fprintf('cancelled\n');
	return
end
% assign name to OptoFileName object to assist name generation
tmpf = OptoFileName(fullfile(plxpath, plxfile));
% build nexinfo name skeleton to search for
nexinfofile = [tmpf.fileWithoutOther '*_nexinfo.mat'];
% get neinfofile and path
[nexinfofile, nexinfopath] = uigetfile(fullfile(plxpath, nexinfofile), ...
													'Select nexinfo.mat file');
