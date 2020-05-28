clc
clear all
warning off

% collect names
current = pwd;
cd(directory)
Type = {to add};
Reference = {'All'};
outfile = to add;

list = dir;
length_names = length(list);
count = 1;
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
        Data.Reference(x,1) = Reference;
        Data.Area_ha(x,1) = data.AREA./10000;
        Data.Mean(x,1) = data.MEAN;
        Data.Maximum(x,1) = data.MAX;
        Data.Sum(x,1) = data.SUM;
        % Uncomment what is needed
        %Data.TransferValue(x,1) = (data.MEAN./100).* Data.Area_ha(x,1); % Aquaduct and Growth days as original 1-km rasters & LPJ resampled to 1-km2 rasters with values as fractions and multiplied by 1M.
       % Data.TransferValue(x,1) = (data.MEAN.*10).*Data.Area_ha(x,1);% InVest
       % Data.TransferValue(x,1) = data.MEAN.*Data.Area_ha(x,1); %$Benefit transfer as original 1-ha rasters
        %Data.TransferValue(x,1) = data.MAX; %G2G and WaterWorld & Luci
       %  Data.TransferValue(x,1) = round((data.MEAN.*((10000/5625)).*Data.Area_ha(x,1))); %Population size form 75 x 75 m rasters (5625 m2), first translated to per hecatre
       Data.TransferValue(x,1) = data.MEAN; %All comparison layers EXCEPT Pop Density and cells
       %Data.TransferValue(x,1) = data.MEAN./1600; %Cells comparison
       %Data.TransferValue(x,1) = (data.SUM./56.25)./(data.AREA/1000000); % comparison pop_density done from the 75m grid map on 10m resolution
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
        Data.Reference(x,1) = Reference;
        Data.Area_ha(x,1) = NaN;
         Data.Mean(x,1) = NaN;
        Data.Maximum(x,1) = NaN;
        Data.Sum(x,1) = NaN;
        Data.TransferValue(x,1) = NaN;
    end  
    clear data
end
cd(current);
save(outfile,'Data')