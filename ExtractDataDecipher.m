clc
clear all
warning off
load('validation_catchments.mat')
data_UK = dataset(UK_catchments,'Varnames','Catchment');
data_Wales = dataset(Wales_catchments,'Varnames','Catchment');
data_Wales.flow(:,1) = NaN;
data_UK.flow(:,1) = NaN;

% collect names
list = dir;
length_names = length(list);
count = 1;
display('Listing files')
for x = 1:1:length_names
    name_list_temp = cellstr(list(x,1).name);
    if strcmp(name_list_temp, '.') ~= 1 && strcmp(name_list_temp, '..') ~= 1
        b = char(name_list_temp);
        c = {b(length(b)-2:length(b))};
        TF = strcmp(c,'csv');
        if TF == 1
            list_of_stations(count,1) =  name_list_temp;
            count = count + 1;
        end
    end
end
clear x

total_stations = length(list_of_stations);
bad_count = 1;
display('Extracting data')
for x = 1:1:total_stations
    fid = fopen(char(list_of_stations(x)),'r');
   
    % Write away the first line
    InputText=textscan(fid,'%s',1,'delimiter','\n;');
 
    % collect the data
    InputText=textscan(fid,'%s','delimiter','\n;');
    leni = length(InputText{1,1});
    counts = 1;
    contains = 0;
    for i = 2:1:leni
        Inputs = InputText{1,1}{i,1};
        count = 1;
        for y = 1:1:length(Inputs)
            tester = Inputs(y);
            if strcmp(tester,',') == 1
                lister(count) = y;
                count = count + 1;
            end
             if strcmp(tester,'"') == 1
                 contains = 1;
             end
        end
        clear count
        % StationID
        if contains == 1
            if i == 2
                Catchment(x,1) = str2double(Inputs(2:(lister(1)-2)));
            end
            year1 = str2double(Inputs((lister(2)+8):(lister(2)+11)));
            year2 = str2double(Inputs((lister(2)+2):(lister(2)+5)));
        else
            if i == 2
                Catchment(x,1) = str2double(Inputs(1:(lister(1)-1)));
            end
            year1 = str2double(Inputs((lister(1)+7):(lister(1)+10)));
            year2 = str2double(Inputs((lister(1)+1):(lister(1)+4)));
        end
        if isnan(year1)== 0
            year = year1;
        elseif isnan(year2) == 0
            year = year2;
        end
        clear year1
        clear year2
        if year <=2015 && year >= 1996
            array(counts,1) = year; %#ok<*SAGROW>
            maxlength = length(lister);
           for t = 2:1:(maxlength-1) % 100
               array(counts,t) = str2double(Inputs((lister(t)+1):(lister(t+1)-1)));
           end
           %array(counts,t+1) =
           %str2double(Inputs(lister(maxlength+1):length(Inputs))); % last
           %one; NOT NEEDED SINCE ENDING ON "," (different from future
           %flows)
           array(counts,t+1)= (nanmean(array(counts,2:(t)))).*86400; 
           % mean Total flow in m3 per day extrapolated per day from m3 per sec with seconds per day;
           counts = counts + 1;
           divider = t-1;
           clear t
        end 
    end
    clear i
    years = unique(array(:,1));
    
    for i = 1:1:length(years)
       array_year(i,1) = sum(array((array(:,1)==years(i)),(divider+2)));
    end
    clear array
    Catchment(x,2) = nanmean(array_year);
    fclose all
    clc
    display(x)
    display(Catchment(x,1))
    display(Catchment(x,2))
    
    where = find(data_UK.Catchment == (Catchment(x,1)));
    if isempty(where)~= 1
    data_UK.flow(where,1) = Catchment(x,2);
    end
    where = find(data_Wales.Catchment == (Catchment(x,1)));
    if isempty(where)~= 1
    data_Wales.flow(where,1) = Catchment(x,2);
    end    
end   
save('Catchments_Decipher','Catchment'); 
save('Decipher_UK', 'data_UK');
save('Decipher_Wales', 'data_Wales')
