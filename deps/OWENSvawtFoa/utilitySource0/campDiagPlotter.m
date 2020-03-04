function campDiagPlotter(resultsFileName,numModesToPlot,numPerRevLines,minRPMplot,maxRPMplot)
%campDiagPlotter plots a Campbell diagram of selected modal frequencies
% **********************************************************************
% *                   Part of the SNL OWENS Toolkit                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   campDiagPlotter(resultsFileName,numModesToPlot,numPerRevLines,...
%                   minRPMplot,maxRPMplot)
%                    
%   This function plots a Campbell diagram of selected modal frequencies
%   vs. rotor speed with per-rev excitation lines (fan plot)
%
%      input:
%      resultsFileName    = string containing results .mat filename
%      numModesToPlot     = number of modes to plot in campbell diagram
%      numPerRevLines     = number of per-rev excitation lines to plot
%      minRPMplot         = minimum rotor speed (RPM) to begin per-rev lines
%      maxRPMPlot         = maximum rotor speed (RPM) to end per-rev lines
% 
%      output:            (NONE)

    load(resultsFileName); %load results file
    close all;             %close all figure windows
    [r,c]=size(freq);

    numModesToPlot = 2*numModesToPlot; %there are 2 of every mode due to state space representation
    if(numModesToPlot > c)  %if numModesToPlot > available number of modes, this adjusts numModesToPlot
        numModesToPlot = c;
        disp('Number of modes requested to plot is greater than available modes in results file.');
    end
    
    colorstr = {'b','g','m','c','y','b','g','m','c','y'};
    for i=2:2:numModesToPlot
           plot(rotorSpeedArray*60,freq(:,i),colorstr{i/2},'LineWidth',2); %plot mode i at various rotor speeds
           hold on;
    end

    %plot per rev lines
    for i=1:numPerRevLines
    linex=[minRPMplot maxRPMplot];
    liney=[minRPMplot maxRPMplot].*i./60.0;
    plot(linex,liney,'--k');
    pstr=[num2str(i),'P'];
    text(0.95*linex(2),liney(2)+.05+(i-1)*.01,pstr,'FontSize',12)
    end
    grid on;
    xlabel('Rotor Speed (RPM)','FontSize',15);
    ylabel('Frequency (Hz)','FontSize',15);
    set(gca,'FontSize',12);
end