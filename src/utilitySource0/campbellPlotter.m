load freqRotorTimo
close all;
[r,c]=size(freq);

numModesToPlot = 12;
numModesToPlot = 2*numModesToPlot;
if(numModesToPlot > c)
	numModesToPlot = c;
end
for i=2:2:numModesToPlot
       plot(rotorSpeedArray*60,freq(:,i),'k','LineWidth',2);
       hold on;
       
       if(i==2)
          %plot flatwise experimental data

            flatwise = importData('snl34mFlatwiseExperimental.csv');
            flatdata = flatwise.data;
            plot(flatdata(:,1),flatdata(:,2),'*r'); 
            
            leadlag =  importData('snl34mLeadLagExperimental.csv');
            leadlagdata = leadlag.data;
            plot(leadlagdata(:,1),leadlagdata(:,2),'*b'); 
       end
end

%plot per rev lines
numPerRevLines = 5;

for i=1:numPerRevLines
linex=[10 60];
liney=[i/6 i];
plot(linex,liney,'--k');
pstr=[num2str(i),'P'];
text(0.95*linex(2),liney(2)+.05+(i-1)*.01,pstr,'FontSize',12)
end
grid on;
xlabel('Rotor Speed (RPM)','FontSize',15);
ylabel('Frequency (Hz)','FontSize',15);
set(gca,'FontSize',12);



h=legend('OWENS','Flatwise Gauge','Lead-Lag Gauge');
set(h,'Location','northwest');