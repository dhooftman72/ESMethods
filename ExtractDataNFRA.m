clc
clear all
warning off

% collect names
list = dir;
length_names = length(list);
count = 1;
for x = 1:1:length_names
    name_list_temp = cellstr(list(x,1).name)
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

total_stations = length(list_of_stations);
bad_count = 1;
for x = 1:1:total_stations
    fid = fopen(char(list_of_stations(x)),'r');
    if x/10 == ceil(x/10)
    display(x)
    end
    % collect the details
    InputText=textscan(fid,'%s',50,'delimiter','\n;');
    % StationID
     a = char(InputText{1,1}{4,1});
     if x == 1
     Catchment_UK = dataset(cellstr(a(12:length(a))),'Varnames',char('StationsID'));
     else
          Catchment_UK.StationsID(x,1) = cellstr(a(12:length(a)));
     end
    % Station_Name
    a = char(InputText{1,1}{5,1});
    Catchment_UK.StationsName(x,1) = cellstr(a(14:length(a)));
    
    
    for z = 1:1:50
        a = (InputText{1,1}{z,1});
        len = min(length(a),9);
        to_comp = (cellstr(a(1:len)));
        test = strcmp('data,last',to_comp);
        if test == 1
            max_heading = z;
        end
    end
    clear z
    % Start date
    a = char(InputText{1,1}{(max_heading-1),1});
    Catchment_UK.StartDate(x,1) = (cellstr(a(12:length(a))));
    % End date
    a = char(InputText{1,1}{max_heading,1});
    Catchment_UK.EndDate(x,1) = (cellstr(a(11:length(a))));
    last_year = str2double(cellstr(a(11:14)));
    fclose all;
    fid = fopen(char(list_of_stations(x)),'r');
    InputText=textscan(fid,'%s',max_heading+1,'delimiter','\n;');
    
    % collect the data
    InputText=textscan(fid,'%s','delimiter','\n;');
    Input = (InputText{1,1});
    leni = length(Input);
    for i = leni:-1:1
         if isempty(Input{i,1}) == 1
             Input(i,:) = [];
         end
    end
    leni = length(Input);
    counts = 1;
    contains = 0;
    for i = 1:1:leni
        Inputs = Input{i,1};
        for y = 1:1:length(Inputs)
            tester = Inputs(y);
            if strcmp(tester,'-') == 1
                contains = 1;
            end
        end
        if contains == 0
            year = str2double(Inputs(7:10));
        else
            year = str2double(Inputs(1:4));
        end
        if year <=2015 && year >= 1996
            len = length(Inputs);
            array(counts,1) = year;
            array(counts,2) = (str2double(Inputs(12:len))).*86400; % mean Total flow in m3 per day extrapolated per day from m3 per sec with seconds per day;
            counts = counts + 1;
        end
    end
    clear i
   if exist('array') == 1; 
        years = unique(array(:,1));
        for i = 1:1:length(years)
            array_year(i,1) = sum(array((array(:,1)==years(i)),2));
        end
    else
        array_year = NaN;
    end
    clear array
    b = char(list_of_stations(x));
    for i = 1:1:length(b)
        check(i,1) =  (strcmp(b(i), '_'));
        end
       ends = (find(check(:,1) == 1))-1;
    Catchment_UK.Mean20y(x,1) = nanmean(array_year);
    Catchment_UK.Std20y(x,1) = nanstd(array_year); 
    Catchment_UK.Check(x,1) =  str2double(b(1:ends));
    clear array_year
    fclose all
    
    clc
    display(x)
    display(Catchment_UK.Check(x,1))
    display(Catchment_UK.Mean20y(x,1))
end
save('Catchments_flow','Catchment_UK');  