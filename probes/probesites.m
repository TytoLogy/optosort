
%% specifications for E-16 probe from Cambridge NeuroTech
% Store specs in a struct, prb
% index to sites (probably unnecessary)
prb.index = 15:-1:0;
% # of sites on probe
prb.nsites = 16;
% xy coordinates of sites, with x location in 1st column, y location
% (increasing depth as value increases) in 2nd column. units are arbitrary,
% but could be um...
% this is as specd in file CNT_E1_Klusta.txt provided by Tahl Holtzman of
% Cambridge NeuroTech 
prb.xy_klusta = [ ...
-20, 160; ...
21, 170; ...
-22, 180; ...
23, 190; ...
-24, 200; ...
25, 210; ...
-26, 220; ...
27, 230; ...
-28, 240; ...
29, 250; ...
-30, 260; ...
31, 270; ...
-32, 280; ...
33, 290; ...
-34, 300; ...
35, 310; ...
];

%% plot the locations as 'o's with site index text adjacent
plot(prb.xy_klusta(:, 1), prb.xy_klusta(:, 2), '.')
text(prb.xy_klusta(:, 1)*(1 + 0.05), prb.xy_klusta(:, 2), num2cell(prb.index))
ylim([150 320])
xlabel('x')
ylabel('y (depth)')
title('raw locations')

%% probe site ID from probe layout pdf: 
%    CambridgeNeuroTech_Acute16ChProbe_Pinouts.pdf
% data also in E1_ProbeLayout.xls
%
% *** these sites are listed from proximal (top of probe) to
% *** distal (tip of probe)
prb.id = [ ...
	9; ...
	11; ...
	13; ...
	10; ...
	12; ...
	15; ...
	14; ...
	16; ...
	8; ...
	6; ...
	4; ...
	2; ...
	1; ...
	3; ...
	5; ...
	7	];

% print out in mock layout
for n = 1:prb.nsites
	if isodd(n)
		fprintf('   %d\n', prb.id(n));
	else
		fprintf('      %d\n', prb.id(n));
	end
end


%% define probe site (column 1) to AD channel (column 2) map:

prb.id_to_channel_map = [ ...
	6,	15; ...
	2,	13; ...
	3,	11; ...
	7,	9; ...
	16,	7; ...
	15,	5; ...
	10,	3; ...
	11,	1; ...
	8,	16; ...
	4,	14; ...
	1,	12; ...
	5,	10; ...
	14,	8; ...
	12,	6; ...
	13,	4; ...
	9,	2 ];

%% for each probe_id site, find the corresponding A/D channel
% so, loop through the sites as listed in prb.id...
prb.channel = zeros(prb.nsites, 1);
for n = 1:prb.nsites
	% find location of this id in prb.site_to_channel_map first column
	indx = find(prb.id(n) == prb.id_to_channel_map(:, 1));
	% matching AD channel will be in second column
	prb.channel(n) = prb.id_to_channel_map(indx, 2);
	fprintf('site %d --> chan %d\n', prb.id(n), ...
						prb.channel(n));
end




%% need to do some tweaking to get sites/channels to line up and 
% match physical layout of probe
% what we need: 
% A/D channels mapped to x, y locations of recording probe sites.

% right now locations of sites in prb.xy_klusta are listed from 
% bottom to top, so flip up down
xyud = flipud(prb.xy_klusta);

% then need to reverse the xy coords in pairs in order to get left-right
% position to match locations on probe
xylr = zeros(size(xyud));
for n = 1:2:prb.nsites
	xylr(n, 1) = xyud(n+1, 1);
	xylr(n+1, 1) = xyud(n, 1);
	xylr(n, 2) = xyud(n, 2);
	xylr(n+1, 2) = xyud(n+1, 2);
end

% check using plot xy coordinates and site id
figure(2)
plot(xylr(:, 1), xylr(:, 2), '.r')
ctxt = cell(length(prb.id), 1);
for n = 1:length(prb.id)
    ctxt{n} = sprintf('id %d', prb.id(n));
    if xylr(n, 1) < 0
        text(xylr(n, 1)*(1 - 0.05), xylr(n, 2), ctxt{n}, 'Color', 'r')
    else
        text(xylr(n, 1)*(1 + 0.05), xylr(n, 2), ctxt{n}, 'Color', 'r')
    end
end
xlim([-40 45])
ylim([150 320])
xlabel('x')
ylabel('y (depth)')
title('tweaked locations');



% Store "tweaked" physical xy locations
prb.xy = xylr;


% Now, need to create map A/D channel to physical location
prb.channel_to_xy_map = [prb.channel prb.xy];

% Create table to hold relevant values
prb.probetable = table(prb.id, prb.channel, prb.xy, 'VariableNames', {'SiteID', 'ADchannel', 'XY'});

% Write to mat file
save('CNT_E1_probestruct.mat', 'prb', '-MAT');

%% Write channel and xy locations in prb format
% This is essentially writing python code to define the geometry of the probe. 


% example:
%{
# header info

total_nb_channels = 30
radius            = 200
channel_groups    = {}

channel_groups = {
    # Shank index.
    1:
        {   
            'channels': list(range(30)),
            'geometry': {
                0: (0, -1050),
                1: (22, -1050),
                2: (0, -1072),
                3: (22, -1072),
                4: (0, -1094),
                5: (22, -1094),
                6: (0, -1116),
                7: (22, -1116),
                8: (0, -1138),
                9: (22, -1138),
                10: (0, -1160),
                11: (22, -1160),
                12: (0, -1182),
                13: (22, -1182),
                14: (0, -1204),
                15: (22, -1204),
                16: (0, -1226),
                17: (22, -1226),
                18: (0, -1248),
                19: (22, -1248),
                20: (0, -1270),
                21: (22, -1270),
                22: (0, -1292),
                23: (22, -1292),
                24: (0, -1314),
                25: (22, -1314),
                26: (0, -1336),
                27: (22, -1336),
                28: (0, -1358),
                29: (22, -1358),
            }
    }
}
%}

% probe file filename
prbfile = 'CNT_E1.prb';
% prbfile = '';

if isempty(prbfile)
    fp = 1;
else
    % open file
    fp = fopen(prbfile, 'w');
end

% define tb (tab) in spaces (3)
tb = '   ';
% write header
fprintf(fp, '# SpyKING CIRCUS .prb configuration file\n');
fprintf(fp, '# for Cambridge NeuroTech E1 16 channel probe\n');
fprintf(fp, '# Date: %s\n', date);

% write # of channels
fprintf(fp, 'total_nb_channels = %d\n', prb.nsites);
% write radius (default)
fprintf(fp, 'radius = %d\n', 200);


fprintf(fp, 'channel_groups = {\n');

fprintf(fp, '%s1: {\n', tb);
fprintf(fp, '%s%s''channels'': list(range(%d)),\n', tb, tb, prb.nsites);
fprintf(fp, '%s%s''geometry'': {\n', tb, tb);

% need to get data sorted by AD channel
srtprobe = sortrows(prb.probetable, 'ADchannel');
for n = 1:prb.nsites
    % need to offset AD channel by 1
    fprintf(fp, '%s%s%s%d: [%d, %d],\n', tb, tb, tb, ...
                        srtprobe.ADchannel(n) - 1, srtprobe.XY(n, 1), srtprobe.XY(n, 2));
end
fprintf(fp, '%s%s}\n', tb, tb);
%fprintf(fp, '%s%s%s''graph'' : [],\n', tb, tb, tb);
fprintf(fp, '%s}\n', tb);
fprintf(fp, '}\n');

% close the file
if fp ~= 1
    fclose(fp);
end











%{



% specifications for E-16 probe from Cambridge NeuroTech
prb.index = 15:-1:0;
prb.nsites = 16;

prb.xy = [ ...
				-20, 160; ...
				21, 170; ...
				-22, 180; ...
				23, 190; ...
				-24, 200; ...
				25, 210; ...
				-26, 220; ...
				27, 230; ...
				-28, 240; ...
				29, 250; ...
				-30, 260; ...
				31, 270; ...
				-32, 280; ...
				33, 290; ...
				-34, 300; ...
				35, 310; ...
			];

% plot the locations as 'o's with site index text adjacent
% this is as specd in file CNT_E1_Klusta.txt
plot(prb.xy(:, 1), prb.xy(:, 2), '.')
text(prb.xy(:, 1)*(1 + 0.05), prb.xy(:, 2), num2cell(prb.index))
ylim([150 320])
%% need to do some tweaking


% probe site ID from probe layout pdf: 
%    CambridgeNeuroTech_Acute16ChProbe_Pinouts.pdf
% data also in E1_ProbeLayout.xls
%{
 *** these sites are listed from proximal (top of probe) to
 *** distal (tip of probe)
  9
      11
   13
      10
   12
      15
   14
      16
   8
      6
   4
      2
   1
      3
   5
      7
    
%}
prb.id = [ ...
	9; ...
	11; ...
	13; ...
	10; ...
	12; ...
	15; ...
	14; ...
	16; ...
	8; ...
	6; ...
	4; ...
	2; ...
	1; ...
	3; ...
	5; ...
	7	];

% print out in mock layout
for n = 1:prb.nsites
	if isodd(n)
		fprintf('   %d\n', prb.id(n));
	else
		fprintf('      %d\n', prb.id(n));
	end
end

% map from probe site ID (1st col) to Omnetics connector (2nd col) map
prb.site_to_omnetics_map = [ ...
	9, 
	

% define probe site (column 1) to AD channel (column 2) map:
prb.site_to_channel_map = [ ...
	6,	1; ...
	2,	3; ...
	3,	5; ...
	7,	7; ...
	16,	9; ...
	15,	11; ...
	10,	13; ...
	11,	15; ...
	8,	2; ...
	4,	4; ...
	1,	6; ...
	5,	8; ...
	14,	10; ...
	12,	12; ...
	13,	14; ...
	9,	16 ];



% for each probe_id site, find the corresponding A/D channel
% so, loop through the sites as listed in prb.id
prb.channel = zeros(prb.nsites, 1);
for n = 1:prb.nsites
	% find location of this id in prb.site_to_channel_map first column
	indx = find(prb.id(n) == prb.site_to_channel_map(:, 1));
	% matching AD channel will be in second column
	prb.channel(n) = prb.site_to_channel_map(indx, 2);
	fprintf('site %d --> chan %d\n', prb.id(n), ...
						prb.channel(n));
end

% what we need: 
% A/D channels mapped to x, y locations of recording probe sites.
% right now locations are listed from bottom to top, so flip up down
xyud = flipud(prb.xy)
% then need to reverse the xy coords in pairs
xylr = zeros(size(xyud));
for n = 1:2:prb.nsites
	xylr(n, 1) = xyud(n+1, 1);
	xylr(n+1, 1) = xyud(n, 1);
	xylr(n, 2) = xyud(n, 2);
	xylr(n+1, 2) = xyud(n+1, 2);
end
xylr

% first: determine channel that matches each site in site.probe_id
figure(2)
plot(xyud(:, 1), xyud(:, 2), '.b')
hold on
plot(xylr(:, 1), xylr(:, 2), '.r')
hold off
ctxt = cell(length(prb.id), 1);
for n = 1:length(prb.id)
	ctxt{n} = sprintf('c%d', prb.id(n));
end
text(xylr(:, 1)*(1 + 0.05), xylr(:, 2), ctxt, 'Color', 'r')
ylim([150 320])

%{
actual coords:
Coordinates	
x	y
-34	310
35	300
-32	290
33	280
-30	270
31	260
-28	250
29	240
-26	230
27	220
-24	210
25	200
-22	190
23	180
-20	170
21	160
%}




%}