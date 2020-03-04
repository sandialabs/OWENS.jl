function plotDispTimeHistory(time,uHist,nodeNum,localDofNum,ylabstr)

numDofPerNode = 6;
globalDOFNum = (nodeNum-1)*numDofPerNode + localDofNum;

plot(time,uHist(globalDOFNum,:));
grid on;
xlabel('time (s)');
ylabel(ylabstr);

end