function RegressAnalysesCarbon
if matlabpool('size') ~= 0
    evalc('matlabpool close');
end
% Do once then comment since loaded below
% Ensembles = ReadEnsembles;
% Variables = ReadVariables;
% [Ensembles, Variables] = CleanDataSetCarbon(Ensembles,Variables);
% save('TheDataperKM','Ensembles','Variables')
%% Regression analyses
clear all
EnsemblesNames = [{'CrossAIC'},{'CrossCross'},{'CrossDev'},...
    {'CrossLeast'},{'CrossRho'},{'UniqueDown'},{'UniqueUp'},...
    {'GridSize'},{'Mean'},{'Median'},{'PCA'},{'Regress2Median'},...
    {'SEMamongEnsembles'},{'STD_among_Models'},{'CVamongEns'},{'CVamongModels'}];
VariableNames = [{'RainfallAnnual'},{'PETAnnual'},{'PopSizeLog'},{'AgriculuturalPerc'},{'GrasslandPerc'},...
    {'PeatPerc'},{'UrbanPerc'},{'WoodlandPerc'},{'TemperatureMeanAnnual'},{'TemperatureMin'},...
    {'TemperatureSeasonality'},{'TemperatureRange'},{'RainfallSeasonality'},{'DEMMean'},{'SlopesMean'}];
VariList = 1:length(VariableNames);
run_max = 10000;
matlabpool open 'Full'
for ensem = 15%length(EnsemblesNames)
    load('TheDataperKM.mat')
    clear TargetEnsemble
    TargetEnsemble = eval( ['Ensembles.',char(EnsemblesNames(ensem))]);
    for vari = VariList
        str = sprintf('      Running variable %s in Ensemble %s',char(VariableNames(vari)),char(EnsemblesNames(ensem)));
        clc
        disp(str)
        obssnames = [{'Longitude'},{'Lattitude'},{'Nr_models'},{'LongxLat'},VariableNames(vari),{'Error'}];
        clear TargetVariable
        TargetVariable =  eval( ['Variables.',char(VariableNames(vari))]);
        if min(TargetVariable)<0
            TargetVariable = TargetVariable+min(TargetVariable);
        end
        TarVar = normaVar(TargetVariable);
        LongitudeNorm = normaVar(Variables.Longitude);
        LattitudeNorm = normaVar(Variables.Lattitude);
        parfor i = 1:run_max
            draw = randperm(length(Variables.Longitude),519);
            Y = TargetEnsemble(draw); %#ok<PFBNS>
            Longitude = LongitudeNorm(draw); %#ok<PFBNS>
            Lattitude = LattitudeNorm(draw); %#ok<PFBNS>
            NrofModels = Variables.NrModels(draw);
            % correcting factors individual R2
            dsInd = dataset(Y,Longitude);
            mdltmp = LinearModel.fit(dsInd,'Y ~ Longitude');
            R_long = min(max(mdltmp.Rsquared.Adjusted,0),1);
            %clear mdltmp dsInd 
            dsInd = dataset(Y,Lattitude);
            mdltmp = LinearModel.fit(dsInd,'Y ~ Lattitude');
            R_lat = min(max(mdltmp.Rsquared.Adjusted,0),1);
            %clear mdltmp dsInd 
            dsInd = dataset(Y,NrofModels);
            mdltmp = LinearModel.fit(dsInd,'Y ~ NrofModels');
            R_Models = min(max(mdltmp.Rsquared.Adjusted,0),1);
            %clear mdltmp dsInd   
            
            % correcting models
            ds = dataset(Y,Longitude,Lattitude,NrofModels);
            mdl = LinearModel.fit(ds,'Y ~ NrofModels + Lattitude+ Longitude + Longitude:Lattitude');
            table = anova(mdl,'component');
            
            % Determinant test
            Yresid = mdl.Residuals.Raw;
            XTarget = TargetVariable(draw); %#ok<PFBNS>
            ds2 = dataset(Yresid,XTarget);
            mdl2 = LinearModel.fit(ds2,'Yresid ~ XTarget');
            tbl2 = anova(mdl2,'component');
            
            SS(:,i)= [double(table(1,1));double(table(2,1));double(table(3,1));double(table(4,1));double(tbl2(1,1));double(tbl2(2,1))]; %#ok<*SAGROW>
            DF(:,i)=[double(table(1,2));double(table(2,2));double(table(3,2));double(table(4,2));double(tbl2(1,2));double(tbl2(2,2))];
            F(:,i)=[double(table(1,4));double(table(2,4));double(table(3,4));double(table(4,4));double(tbl2(1,4));double(tbl2(2,4))];
            P(:,i)=[double(table(1,5));double(table(2,5));double(table(3,5));double(table(4,5));double(tbl2(1,5));double(tbl2(2,5))];
            R(:,i) = [R_long;R_lat;R_Models;0;(min(max(mdl2.Rsquared.Adjusted,0),1));0];
            Direction(:,i)= [mdl.Coefficients.Estimate(2,1);mdl.Coefficients.Estimate(3,1);mdl.Coefficients.Estimate(4,1);...
                mdl.Coefficients.Estimate(5,1);mdl2.Coefficients.Estimate(2,1);0];
        end
        if vari == VariList(1)
            Stats = dataset(nanmean(SS,2),'ObsNames',obssnames','VarNames',{'SumofSquares'});
            Stats.DF = nanmean(DF,2);
            Stats.MeanSquares = nanmean((SS./DF),2);
            Stats.Fvalue = nanmean(F,2);
            Stats.Pvalue = nanmean(P,2);
            Stats.R = nanmean(R,2);
            Stats.Direction = nanmean(Direction,2);
        else
            New = dataset(nanmean(SS(5,:),2),'ObsNames',char(VariableNames(vari)),'Varnames',{'SumofSquares'});
            New.DF = nanmean(DF(5,:),2);
            New.MeanSquares = nanmean((SS(5,:)./DF(5,:)),2);
            New.Fvalue = nanmean(F(5,:),2);
            New.Pvalue = nanmean(P(5,:),2);
            New.R =  nanmean(R(5,:),2);
            New.Direction = nanmean(Direction(5,:),2);
            Stats(char(VariableNames(vari)),:) = New;
        end
        clear SS DF F P R  TargetVariable New
    end
    Results.(genvarname(char(EnsemblesNames(ensem)))) = Stats;
    clearvars -except Results EnsemblesNames VariableNames VariList run_max
end
save('ResultsCarbonNew','Results','EnsemblesNames', 'VariableNames', 'VariList');
matlabpool close
end

function Ensembles = ReadEnsembles
% Read all Variables
% Read list
cd Ensemble layers
list = dir;
length_names = length(list);
count = 1;
for x = 1:1:length_names
    name_list_temp = cellstr(list(x,1).name);
    if strcmp(name_list_temp, '.') ~= 1 && strcmp(name_list_temp, '..') ~= 1
        b = char(name_list_temp);
        c = {b(length(b)-2:length(b))};
        TF = strcmp(c,'asc');
        if TF == 1
            list_of_ensembles(count,1) =  name_list_temp;
            count = count + 1;
        end
    end
end

for i = 1:(count-1)   
    file = char(list_of_ensembles(i));
    for x = 1:1:length(file)
        check(x) =  (strcmp(file(x), '.'));
    end
    ends = (find(check == 1))-1;
    VarName = file(1:ends);
    clear check x
    if i == 1
        [data_block,ncols,nrows] = cut_top_off_ascii(file);
        Var = reshape(data_block,ncols*nrows,1);
        Var(Var == -9999) = NaN;
        Ensembles = dataset(Var,'VarNames',char(VarName));
    else
        [data_block,~,~] = cut_top_off_ascii(file);
         Var = reshape(data_block,ncols*nrows,1);
          Var(Var == -9999) = NaN;
          Ensembles.(genvarname(char(VarName)))= Var;
    end
end
cd ..
end


function  Variables = ReadVariables
% Read list
list = dir;
length_names = length(list);
count = 1;
for x = 1:1:length_names
    name_list_temp = cellstr(list(x,1).name);
    if strcmp(name_list_temp, '.') ~= 1 && strcmp(name_list_temp, '..') ~= 1
        b = char(name_list_temp);
        c = {b(length(b)-2:length(b))};
        TF = strcmp(c,'asc');
        if TF == 1
            list_of_variables(count,1) =  name_list_temp;
            count = count + 1;
        end
    end
end

for i = 1:(count-1)   
    file = char(list_of_variables(i));
    for x = 1:1:length(file)
        check(x) =  (strcmp(file(x), '.'));
    end
    ends = (find(check == 1))-1;
    VarName = file(1:ends);
    clear check x
    if i == 1
        [data_block,ncols,nrows] = cut_top_off_ascii(file);
        Var = reshape(data_block,ncols*nrows,1);
        Var(Var == -9999) = NaN;
        Variables = dataset(Var,'VarNames',char(VarName));
    else
        [data_block,~,~] = cut_top_off_ascii(file);
         Var = reshape(data_block,ncols*nrows,1);
          Var(Var == -9999) = NaN;
          Variables.(genvarname(char(VarName)))= Var;
    end
end
end

function [data_block,ncols,nrows] = cut_top_off_ascii(file)
number_of_heading_lines = 6;
fid = fopen(file,'r');
InputText=textscan(fid,'%s',number_of_heading_lines,'delimiter','\n');
for head = 1:1:number_of_heading_lines
    temp = InputText{1,1}{head,1};
    length_head = length(temp);
    temp_text = temp(1,15:length_head);
    Intro(head) = str2num(temp_text);
    clear temp
    clear length_head
    clear temp_text
end
ncols = Intro(1);
nrows = Intro(2);
xllcorner = Intro(3);
xllcorner = round(xllcorner);
yllcorner = Intro(4);
yllcorner = round(yllcorner);
cellsize =   Intro(5);
NODATA_value =  Intro(6);

Block = 1;   % Initialize block index
sprintf('Block: %s', num2str(Block));                % Display block number
InputText=textscan(fid,'Num SNR=%f'); % Read parameter value
NumCols=InputText{1};
FormatString=repmat('%f',1,NumCols);  % Create format string based on parameter
InputText=textscan(fid,FormatString,'delimiter',','); % Read data block
Data{Block,1}=cell2mat(InputText); %#ok<*SAGROW> % Convert to numerical array from cell
data_array = Data{1,1};
count = 1;
for row = 1:1: nrows
    for col = 1:1:ncols
        data_block(row,col) = data_array(count); %#ok<*AGROW>
        count =count + 1;
    end
end
fclose('all');
end

function [Ensembles, Variables] = CleanDataSetCarbon(Ensembles,Variables)
Ensembles_dbl = double(Ensembles);
Variables_dbl = double(Variables);
test = find(isnan((Ensembles_dbl(:,1))) == 1);
Ensembles_dbl(test,:) = [];
Variables_dbl(test,:) = [];
Ensembles_new = dataset(Ensembles_dbl(:,1),'Varnames','CrossAIC');
Ensembles_new.CrossCross = Ensembles_dbl(:,2);
Ensembles_new.CrossDev = Ensembles_dbl(:,3);
Ensembles_new.CrossLeast = Ensembles_dbl(:,4);
Ensembles_new.CrossRho = Ensembles_dbl(:,5);
Ensembles_new.UniqueDown = Ensembles_dbl(:,6);
Ensembles_new.UniqueUp = Ensembles_dbl(:,13);
Ensembles_new.GridSize = Ensembles_dbl(:,7);
Ensembles_new.Mean = Ensembles_dbl(:,8);
Ensembles_new.Median = Ensembles_dbl(:,9);
Ensembles_new.PCA = Ensembles_dbl(:,10);
Ensembles_new.Regress2Median = Ensembles_dbl(:,11);
Ensembles_new.SEMamongEnsembles = Ensembles_dbl(:,12);
Ensembles_new.Lattitude = Variables_dbl(:,5);
Ensembles_new.Longitude = Variables_dbl(:,6);
Variables_new = dataset(Variables_dbl(:,6),'Varnames','Longitude');
Variables_new.Lattitude = Variables_dbl(:,5);
Variables_new.CarbonZones = Variables_dbl(:,3);
Variables_new.RainfallAnnual = Variables_dbl(:,15);
Variables_new.PETAnnual = Variables_dbl(:,12);
Variables_new.PopSize = Variables_dbl(:,13);
Variables_new.AgriculuturalPerc = Variables_dbl(:,1)./1600;
Variables_new.GrasslandPerc = Variables_dbl(:,4)./1600;
Variables_new.PeatPerc = Variables_dbl(:,11)./1600;
Variables_new.UrbanPerc = Variables_dbl(:,19)./1600;
Variables_new.WoodlandPerc = Variables_dbl(:,20)./1600;
Variables_new.TemperatureMeanAnnual = Variables_dbl(:,7);
Variables_new.TemperatureMin = Variables_dbl(:,9);
Variables_new.TemperatureSeasonality = Variables_dbl(:,18);
Variables_new.TemperatureRange = Variables_dbl(:,16);
Variables_new.RainfallSeasonality = Variables_dbl(:,14);
Variables_new.NrModels = Variables_dbl(:,10);
Variables_new.DEMMean = Variables_dbl(:,8);
Variables_new.SlopesMean = Variables_dbl(:,17);
Ensembles = Ensembles_new;
Variables = Variables_new;
end
