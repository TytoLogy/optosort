function varargout = tracePlot(pH, T, Fs)
   % figure and axes background grey level ( 0 = black, 1 = white)
   bgcol = 1;

   % max range for each trace
   yabsmax = 5;

   [npts, nchan] = size(T);
   
   if isempty(pH)
      % setup plot
      % generate figure, axes
      fH = figure;
      aX = axes;
      % switch focus to figure
      figure(fH);
      % time vector for xdata
      tv = timevec(npts, Fs, 'ms');
      % generate temporary data
      tmpData = zeros(size(T));
      for n = 1:nchan
         tmpData(:, n) = n*(yabsmax) + 2*(2*rand(npts, 1)-1);
      end
      % and plot it
      pH = plot(aX, tv, tmpData, 'k-');
      yticks_yvals = yabsmax*(1:nchan);
      yticks_txt = cell(nchan, 1);
      for n = 1:nchan
         yticks_txt{n} = num2str(n);
      end
      ylim(yabsmax*[0 nchan+1]);
      set(aX, 'YTick', yticks_yvals);
      set(aX, 'YTickLabel', yticks_txt);
      set(aX, 'TickDir', 'out');
      set(aX, 'Box', 'off');
   %    set(fH, 'Position', [792 225 557 800]);		
      xlabel('Time (ms)')
      ylabel('Channel')
      set(aX, 'Color', bgcol*[1 1 1]);
      set(fH, 'Color', bgcol*[1 1 1]);
      set(fH, 'ToolBar', 'none');
      varargout{1} = fH;
      varargout{2} = aX;
      varargout{3} = pH;
   else
	   % assign data to plot
	   for c = 1:nchan
		   % update plot
		   set(pH(c), 'YData', T(:, c) + c*yabsmax);
	   end
	   % update plots
	   refreshdata(ancestor(pH(1), 'figure'));
   end

end