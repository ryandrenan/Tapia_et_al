%% GENERAL MED PC to MATLAB Script for IVSA Data
% 
% version: 7
% date: 12.01.21
% minor update: n/a 

%% Process names of files to run
% point matlab to directory containing desired mpc files
format short g %makes numbers easier format to read on screen
close all %close all existing graphs
clear all %clear all vars

% can change to default finder window to elsewhere (currently Desktop)
batch_dir = '~/Desktop/';
dirinfo_ = uigetdir(batch_dir,'Directory');
dirinfo = dir(dirinfo_);
dirinfo([dirinfo.isdir]) = [];  % remove directories

% remove any non-data files assuming only data files have the word 'Subject'
dirinfo(~contains({dirinfo.name},{'Subject'})) = [];

%% Run through dirinfo list, opens each file and parses it
for w = 1:size(dirinfo,1)
    mpc = dirinfo(w).name
    fid = fopen(mpc,'r'); 
    file = fread(fid, 'uint8=>char')'; %or can fscanf(fileID, '%c')

    % splits by new line
    file = splitlines(file); %turns into cell array
    file((file == "")) = []; %removes empty lines

    % find C Array because that has information we want to process
    C = find(file == "C:"); % finds index where C Array begins
    Cend = file(end); % save temporary copy of the last row of C Array
    
    %for some weird reason, this particular case makes C get flipped. this
    %is when there are exactly 4 values in the C array, where there is one
    %value on the 3: line
    if (C+1 == length(file)-1)
        C = split(strip(file(C+1:end-1))); 
        C = C'; %have to flip the vector back to the right dimension
    else
        C = split(strip(file(C+1:end-1))); %most of the time, this runs. extracts C Array without last row
    end
    
    % now splits last row - data files don't always have 3 columns recorded.
    Cend = split(strip(Cend))';
        
    % expand last row of C to have 3 columns of data (should be 1 x 4 cell)
    if size(Cend,2) < 4
        i = 4-size(Cend,2);
        
        for j = 1:i
            Cadd(1,j) = [NaN]; % makes temporary 1-row var w/ columns to add
        end
               
        Cend = [Cend num2cell(Cadd)];
    end
    
    % reattach last row to C Array data
    %if Cend is nan and if it is also index 5 or greater, cut it off

    C = [C; Cend];

    % now we want to remove first column
    % (it shows array element # recorded by medpc, e.g. "0:, 3:, ...")
    % AND reshape into vector while preserving order (important!)
    C = str2double(C(:,2:end));
    C = reshape(C',1,[])';
    C(isnan(C))=[]; %removes NaN before continuing
    rawresults(w).Carray = C; %so we can extract the needed info from the loop process for analysis below
    clearvars Cadd; %needed when using this code in a loop for batch
end
%end of parsing the mpc file, next section is processing it into a time and
%results vector
%the rawresults struct can be exported here to a file if desired, to have a
%copy of the raw mpc data

%% Process the raw data in time and results vector
w = size(rawresults,2);
for p = 1:w
    ten = rawresults(p).Carray
            
    %use any function to remove NaN values
    tem = ten(any(ten,2),:);

    %length of current column of data
    len = length(tem);
    
    %initialize temporary times and results vector of length of current
    %column
    times = zeros(len,1);
    results = zeros(len,1);

    %make a special case here for when len = 2, the smallest it could be
    %(no responses)
    if len == 2
        times(1) = 1;
        times(2) = 1;
        results(1) = 0.7;
        results(2) = 0.7;
    else
        
    %loop that finds numbers greater than 1 and treats them as time values
    %and finds numbers less than 1 and treats them as results values,
    %placing the large/small numbers into the temporary times or results
    %vector
    for i = 1:(len-1)
        if ((tem(i) > 1) & (tem(i+1) < 1))
            times(i) = floor(tem(i)); %floor cuts off the decimal and rounds down
            results(i) = tem(i+1);
        end
        if ((tem(i) > 1) & (tem(i+1) > 1))
            times(i) = floor(tem(i));
            results(i) = 0.4; %0.4 is assigned when there is a timeout press
        end
    end
    end
    %the times and results vectors generated above will have zeros
    %throughout them since it was initialized as a zero vector, so the
    %following line of code removes those zeros
    a = times(any(times,2),:);

    %the MedPC values have the first time value as the first lever
    %response, and each subsequent time value is relative to the first time
    %value, so the following code makes a new vector with cumulative times
    %from the start to the end of the session
    inttimes = a; %intermediate time vector
    inttimes(1) = a(1); %int time vector, first value is that special first value
    for q = 2:length(a)
        inttimes(q) = a(q) + inttimes(q-1);
    end

    %MedPC time values need to be multiplied by 0.01 sec to be converted to
    %seconds
    newtimes = inttimes * 0.01;
    
    %as above, remove any zeroes from the newresults vector
    newresults = results(any(results,2),:);
    
    %initialize a struct array with 2 fields: times and results, where w,
    %the current column from excel, is the counter for which struct in the
    %array i am in. there should wind up being the same number of structs
    %as there were columns of data imported from excel
    ivsatrials(p).times = newtimes;
    ivsatrials(p).results = newresults;
    x = find(ivsatrials(p).results == 0.2);
    ivsatrials(p).infusiontimes = ivsatrials(p).times(x);
    g = find(ivsatrials(p).results == 0.5);
    ivsatrials(p).inactivetimes = ivsatrials(p).times(g);
    
    inftimescurrent = ivsatrials(p).infusiontimes;
    numinf = numel(inftimescurrent); %determine number of infusions
    inflat = zeros(numinf,1); %initialize infusion latency vector for current trial
    for d = 2:numinf
        inflat(d,1) = inftimescurrent(d) - (20+inftimescurrent(d-1));
    end
    inflatnz = inflat(find(inflat > 0));
    ivsatrials(p).infusionlatency = inflatnz;
    ivsatrials(p).meaninfusionlatency = mean(inflatnz);
    
    resnum = length(newresults); %number of responses in current trial
    
    %initialize a few vars for calculating things
    infusions = 0; 
    TOresponse = 0;
    inactive = 0;
    dose = 0.03; % in mg/kg/inf...might want to have this supplied by the user instead
    
    % extract results from current trial
    for b = 1:resnum
        if newresults(b) == 0.2
            infusions = infusions + 1;
        end
        if newresults(b) == 0.4
            TOresponse = TOresponse + 1;
        end
        if newresults(b) == 0.5
            inactive = inactive + 1;
        end
    end
    
    %calculate and assign to ivsatrials structure for later plotting
    ivsatrials(p).active = infusions + TOresponse;
    ivsatrials(p).infusions = infusions;
    ivsatrials(p).inactive = inactive;
    ivsatrials(p).NICintake = (infusions * dose);
    
    %extract session duration for each session
    tempo = strip(file(26)); % cell 26 of the 'file' variable should contain the line indicating the number of minutes that the session duration is.
    leng = split(strip(tempo)); %this separates the line indicator 'S:' from the number
    ivsatrials(p).sessionduration = round(str2double(leng(2))); %round converts to integer and str2double(leng(2)) converts the string at index 2 to a double
end
%%


%% save ivsatrials to mat file and excel
save ivsatrials
T = struct2table(ivsatrials,'AsArray',true);
filename = 'ivsadata.xlsx';
writetable(T,filename);

%% Identify Session Duration
dura = ivsatrials(1).sessionduration; % uses duration of first session, assumes all session durations are the same

%% Plots of the data 

%% Hisogram of Latency Values
figure
alllatencies = vertcat(ivsatrials.infusionlatency);
hold on
histogram(alllatencies,50);
xlabel('Latency (sec)')
ylabel('Number')
title('Latency to Next Infusion after TO');
hold off
save alllatencies

filename = 'alllatencies.xlsx';
writematrix(alllatencies,filename);

%% Raster plot 
figure
box off
hold on
tiledlayout(2,1) %two plots, one upper, one lower
nexttile
%raster plot of all response types

% for loop for which trial we are plotting
for j = 1:length(ivsatrials)
    
% assign the current trial time and result vector to working variables
t1 = ivsatrials(j).times;
r1 = ivsatrials(j).results;

%the following plots a light yellow box under the raster plot, for contrast
%the box goes from (0,0) to (7200,0) to (7200,ymax) to (0,ymax)
x1 = [0 dura dura 0];
y1 = [0 0 length(ivsatrials) length(ivsatrials)];
if j == 1 %only need to plot it once, before any lines are plotted
    patch(x1,y1,[1 0.98 0.839]) %find hex color, convert to RBG, then scale
end

%plotting for loop, plots a tick for each lever press, color coded based on
%whether it was reinforced, was during the TO, or was on the inactive lever
    for i = 1:length(t1)
        %plots line from point (0,t1(i)) to (1, t1(i)) of color and width
        %the 0.1 is added as an offset so there is space between the rows
        %of lines
        if r1(i) == 0.2 %reinforced response is red
            line([(t1(i)) (t1(i))], [((j-1)+0.1) (j-0.1)], 'Color', 'red', 'LineWidth', 1)
        end
        if r1(i) == 0.4 %during TO is blue
            line([(t1(i)) (t1(i))], [((j-1)+0.1) (j-0.1)], 'Color', 'blue', 'LineWidth', 0.5)
        end
        if r1(i) == 0.5 %inactive lever 
            line([(t1(i)) (t1(i))], [((j-1)+0.1) (j-0.1)], 'Color', 'black', 'LineWidth', 1)    
        end
    end
end

xlim([0 dura])
ylim([0 length(ivsatrials)]) %sets y max to the last trial number
xlabel('Time (min)')
ylabel('Training Session #')
xticks([0:1200:dura])

set(gca,'TickDir','out');
title('All Active/Inactive Responses');
hold off

%raster plot of only infusions
nexttile
hold on

for j = 1:length(ivsatrials)
    
%the following plots a light yellow box under the raster plot, for contrast
%the box goes from (0,0) to (14400,0) to (14400,ymax) to (0,ymax)
x1 = [0 dura dura 0];
y1 = [0 0 length(ivsatrials) length(ivsatrials)];
if j == 1 %only need to plot it once, before any lines are plotted
    patch(x1,y1,[1 0.98 0.839]) %find hex color, convert to RBG, then scale
end

newt1 = ivsatrials(j).infusiontimes; % place infusion times for current trial in a working variable for plotting

%plotting for loop, plots a tick for each infusion
    for i = 1:length(newt1)
        line([(newt1(i)) (newt1(i))], [((j-1)+0.1) (j-0.1)], 'Color', 'red', 'LineWidth', 1)
    end
end
grid on
xlim([0 dura])
ylim([0 length(ivsatrials)]) %sets y max to the last trial number
xticks([0:1200:dura])

xlabel('Time (min)')
ylabel('Training Session #')
set(gca,'TickDir','out');
title('Infusions');
hold off


%% Cumulative response graph for active vs. inactive responses
%only plots it for the most recent trial

last = length(ivsatrials);
newt3 = ivsatrials(last).times;
newr3 = ivsatrials(last).results;

for i = 1:length(newt3) %go thru newt3
    if newr3(i) == 0.5 %if an inactive press, zero it to be deleted
        newt3(i) = 0; 
        newr3(i) = 0;
    end
        
%then go through that resultant vector and remove any zeros to arrive at a
%collapsed times vector representing only the infusion times
    newt4 = newt3(any(newt3,2),:);
    newr4 = newr3(any(newr3,2),:);
end


figure
hold on
%draw cumulative response graph for active presses

if newr4(1) == 0.7
    line ([0 dura],[0 0], 'Color','blue','LineWidth',0.75)
else

%first horizontal
line([0 newt4(1)],[0 0], 'Color','blue','LineWidth',0.75)

%vertical is contingent on press type
if newr4(1) == 0.2
    line([newt4(1) newt4(1)],[0 1], 'Color','red','LineWidth',2)
else
    line([newt4(1) newt4(1)],[0 1], 'Color','blue','LineWidth',0.75)
end

for c = 1:length(newt4)
    if c < length(newt4)
        %horizontal
        line([newt4(c) newt4(c+1)], [c c],'Color','blue','LineWidth',0.75)
        %vertical
        if newr4(c+1) == 0.2
            line([newt4(c+1) newt4(c+1)],[c (c+1)], 'Color','red','LineWidth',2)
        else
            line([newt4(c+1) newt4(c+1)],[c (c+1)], 'Color','blue','LineWidth',0.75)
        end
    end
    if c == length(newt4)
        %horizontal
        line([newt4(c) dura], [c c], 'Color','blue','LineWidth',0.75)
    end
end

%draw cumulative response graph for inactive lever presses to be overlaid
%on the infusion plot from above

newt2 = ivsatrials(p).inactivetimes; %newt2 is the time vector for inactive presses

if length(newt2) ~= 0
    line([0 newt2(1)],[0 0], 'Color','black','LineWidth',0.75)
    line([newt2(1) newt2(1)],[0 1], 'Color','black','LineWidth',0.75)
    for g = 1:length(newt2)
        if g < length(newt2)
            line([newt2(g) newt2(g+1)], [g g],'Color','black','LineWidth',0.75)
            line([newt2(g+1) newt2(g+1)],[g (g+1)], 'Color','black','LineWidth',0.75)
        end
        if g == length(newt2) %if at the last response, plot a finishing line
            line([newt2(g) dura], [g g], 'Color','black','LineWidth',0.75)
        end
    end
end

if length(newt2) == 0
    g = 0;
end

xlim([0 dura])
if c > g
    ylim([0 c+1])
end
if g > c
    ylim([0 g+1])
end
ax = gca;
ax.FontSize = 13;
title('Cumulative Responses, Most Recent Session')
xlabel('Time (minutes)')
ylabel('Cumulative Responses')
xticks([0:1200:dura])


end
grid on
hold off

%% Plot responses by day, infusions, and intake

activepokes = zeros(last,1);
inactivepokes = zeros(last,1);
nicinfusions = zeros(last,1);
intake = zeros(last,1);

for s = 1:last
    activepokes(s) = ivsatrials(s).active;
    inactivepokes(s) = ivsatrials(s).inactive;
    nicinfusions(s) = ivsatrials(s).infusions;
    intake(s) = ivsatrials(s).NICintake;
end

figure

tiledlayout(3,1)
nexttile
hold on
grid on
plot(activepokes, '-o', 'Color', 'red')
plot(inactivepokes, '-o', 'Color', 'black')
title('Active/Inactive Responses')
xlabel('Session #')
ylabel('Responses')
xlim([0 (last+1)])
maxapokes = max(activepokes);
maxipokes = max(inactivepokes);
if maxapokes > maxipokes
    ylim([0 (maxapokes*1.1)])
else
    ylim([0 (maxipokes*1.1)])
end
xticks([1:2:last])

hold off

nexttile
hold on
grid on
plot(nicinfusions, '-o', 'Color', 'blue')
line([0 last], [10 10],'Color', 'red', 'LineStyle', '--', 'LineWidth', 1)
title('Drug Infusions')
xlabel('Session #')
ylabel('Infusions Earned')
xlim([0 (last+1)])
ylim([0 (max(nicinfusions)*1.1)])
xticks([1:2:last])
hold off

nexttile
hold on
grid on
plot(intake, '-o', 'Color', 'magenta')
title('Drug Intake')
xlabel('Session #')
ylabel('Nicotine Intake (mg/kg)')
xlim([0 (last+1)])
ylim([0 (max(intake)*1.1)])
xticks([1:2:last])
hold off



