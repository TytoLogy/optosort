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

% *** these sites are listed from proximal (top of probe) to
% distal (tip of probe) ***
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

% probe site (column 1) to AD channel (column 2) map:
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
text(xylr(:, 1)*(1 + 0.05), xylr(:, 2), num2cell(prb.id))
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