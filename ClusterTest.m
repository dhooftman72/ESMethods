function Parameters = ClusterTest(Parameters)
reply = 'N';
while reply == 'N'
reply = 'Y'; %#ok<NASGU>
str = sprintf('      Uniqueness will be done with %d categories ',Parameters.NrCat);
disp(str)
test = Parameters.UniqueMatrix(Parameters.SetNames,:);
 Y =squareform(pdist(double(test),'spearman'));
Z = linkage(Y);
Parameters.groups = clusterdata(double(test),'distance','spearman','maxclust',Parameters.NrCat);
for i = 1:length(Parameters.SetNames)
    Text =[int2str(Parameters.groups(i)),'_',char(Parameters.SetNames(i))];
   LabelsText(i) = {Text};     %#ok<*AGROW>
   clear Text
end
dendrogram(Z,'labels',LabelsText)
t = timer('ExecutionMode', 'singleShot', 'StartDelay', Parameters.TimeOut, 'TimerFcn', @pressEnter);
start(t);
reply = input('      Is this the right number of categories Y/N: ', 's');
stop(t);
if isempty(reply) == 1
    reply = 'Y';
end
if strcmp('y',reply) == 1
    reply = 'Y';
end
if strcmp('Y',reply) ~= 1 && strcmp('y',reply) ~= 1
    reply = 'N';
end
str = sprintf('      The reply was: %s ',reply);
disp(str)
str = sprintf('      Uniqueness will be done with %d categories ',Parameters.NrCat);
disp(str)
if strcmp('Y',reply) ~= 1
    Parameters.NrCat = str2double(input('Enter new number of categories: ', 's'));
end
close
end
end

function pressEnter(HObj, event) %#ok<INUSD>
  import java.awt.*;
  import java.awt.event.*;
  rob = Robot;
  rob.keyPress(KeyEvent.VK_ENTER)
  rob.keyRelease(KeyEvent.VK_ENTER)
end
