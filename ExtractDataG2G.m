clc
clear all
warning off

% collect names
list = dir;
length_names = length(list);
count = 1;
Type = {'G2G'};
Reference = {'Wales'};
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
clear  name_list_temp

total_stations = length(list_of_stations);
bad_count = 1;
Data = dataset(NaN,'Varnames','Catchment');
for x = 1:1:total_stations
    data = dataset('file',char(list_of_stations(x)),'delimiter',',','ReadObsNames',false);
    if isempty(data) ~= 1
        Data.Catchment(x,1) = data.ID_STRING;
        Data.Type(x,1) = Type;
        Data.Referenceype(x,1) = Reference;
        Data.Area_ha(x,1) = data.AREA./10000;
        Data.Maximum(x,1) = data.MAX;
        Data.Sum(x,1) = data.SUM;
    else
        b = char(list_of_stations(x));
        for i = 1:1:length(b)
            check(i,1) =  (strcmp(b(i), '_'));
        end
        for i = 1:1:length(b)
            check(i,2) =  (strcmp(b(i), '.'));
        end
        start = (find(check(:,1) == 1))+1;
        ends = (find(check(:,2) == 1))-1;
        Data.Catchment(x,1) = str2double(b(start:ends)); 
        Data.Type(x,1) = Type;
        Data.Area_ha(x,1) = NaN;
        Data.Maximum(x,1) = NaN;
        Data.Sum(x,1) = NaN;
    end  
    clear data
end