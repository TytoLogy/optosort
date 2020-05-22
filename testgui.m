function testgui(S)

	% add GUI components
	G = addComponents(S);
	% make figure visible
	G.fig.Visible = 'on';
	
	function G = addComponents(S)
		% Create GUI
		G.fig = figure(	'Units', 'Points', ...
								'Visible', 'off', ...
								'Resize', 'off', ...
								'Tag', 'fig');
		G.plotBtn = uicontrol(	G.fig, ...
									'Position',[10 340 70 30],...
									'String','Plot Sine',...
									'Tag','plotBtn',...
									'Callback',@plotsine);
		
		fnames = cell(S.Info.nFiles, 1);
		for f = 1:S.Info.nFiles
			fnames{f} = S.Info.FileInfo{f}.F.base;
		end
								
		G.fileSelPulldown = uicontrol(G.fig, ...
											'Style', 'popupmenu', ...
											'Position', [10, 300, 100, 30], ...
											'String', fnames, ...
											'Tag', 'fileSelPulldown', ...
											'Value', 1);
											
										
								
								
								
								
		G.ax = axes(	'Parent', G.fig,...
							'Position', [0.20 0.13 0.71 0.75],...
							'Tag','ax');
	end

	function plotsine(hObject, event)
		theta = 0:pi/64:6*pi;
		y = sin(theta);
		plot(G.ax,theta,y);
	end
end



