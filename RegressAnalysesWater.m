function RegressAnalysesWater
% Collate data from Ensemble calculations
% Do once then comment since loaded below
%[VariablesWater,EnsemblesWater] = CollateWaterData;

%% Regression analyses
clear all
if matlabpool('size') ~= 0
    evalc('matlabpool close');
end
EnsemblesNames = [{'Mean'},{'Median'},{'PCA_weighting'},...
    {'Median_among'},{'RegressAmong'},{'MaxentDev'},{'MaxentRho'},...
    {'MaxentLeast'},{'CorCoef'},{'GridSize'},{'UniqueUpWeight'},{'UniqueDownWeight'},...
    {'SEMamongEnsembles'},{'SEM_among_Models'},{'CVamongEns'},{'CVamongModels'}];
VariableNames = [{'PopSizeLog'},{'RainfallAnnual'},{'PETAnnual'},{'RainfallSeasonality'},...
                    {'TemperatureMeanAnnual'},{'TemperatureMin'},{'TemperatureSeasonality'},...
                    {'TemperatureRange'},{'DEMMean'},{'SlopesMean'},... 
                    {'AgriculuturalPerc'},{'GrasslandPerc'},...
                    {'PeatPerc'},{'UrbanPerc'},{'WoodlandPerc'}];
VariList = [1:length(VariableNames)];
for ensem = 13:14%1:length(EnsemblesNames)
    load('TheWaterDataperKM.mat')
    TargetEnsemble = eval( ['EnsemblesWater.',char(EnsemblesNames(ensem))]);
    for vari = VariList
        str = sprintf('      Running variable %s in Ensemble %s',char(VariableNames(vari)),char(EnsemblesNames(ensem)));
        clc
        disp(str)
        obssnames = [{'SpatialAutoCorrect'},VariableNames(vari),{'Error'}];
        TargetVariable =  eval( ['VariablesWater.',char(VariableNames(vari))]);
        if min(TargetVariable)<0
            TargetVariable = TargetVariable+min(TargetVariable);
        end
        if min(TargetVariable)<0
            TargetVariable = TargetVariable+min(TargetVariable);
        end
        TarVar = normaVar(TargetVariable);
       Longitude = VariablesWater.Longitude;
       Lattitude = VariablesWater.Lattitude;
        Y = TargetEnsemble;
        XTarget =  TarVar;
        
        [~,~,~,wij] = Morans(Longitude,Lattitude,250, Y,1); 
            %AutoCorrelation = Autoset(AutoCorrelation,Size,{'OverallAccuracy'},Morans_I,PValue,ZValue,1);
            Auto = zeros(519,1);
            for s = 1:519
                Autot = 0;
                for t = 1:519
                    if s ~= t
                        Autot = Autot + (wij(s,t).*Y(t));
                    end
                end
                Auto(s,1) = Autot./sum(wij(s,:)); %#ok<*AGROW>
            end
            % Determinant test
            
             % correcting models
            ds = dataset(Y,Auto);
            mdl = LinearModel.fit(ds,'Y ~ Auto');
            table = anova(mdl,'component');
            
            Yresid = mdl.Residuals.Raw;
            ds2 = dataset(Yresid,XTarget);
            mdl2 = LinearModel.fit(ds2,'Yresid ~ XTarget');
            tbl2 = anova(mdl2,'component',1);

            SS(:,1)= [double(table(1,1));double(tbl2(1,1));double(tbl2(2,1))]; %#ok<*SAGROW>
            DF(:,1)= [double(table(1,2));double(tbl2(1,2));double(tbl2(2,2))];
            F(:,1)=[double(table(1,4));double(tbl2(1,4));double(tbl2(2,4))];
            P(:,1)=[double(table(1,5));double(tbl2(1,5));double(tbl2(2,5))];
            R(:,1) = [(min(max(mdl.Rsquared.Adjusted,0),1));(min(max(mdl2.Rsquared.Adjusted,0),1));0];
            Direction(:,1)= [mdl.Coefficients.Estimate(2,1);mdl2.Coefficients.Estimate(2,1);0];
            Intercept(:,1) = [mdl.Coefficients.Estimate(1,1);mdl2.Coefficients.Estimate(1,1);0];
        %end
        if vari == VariList(1)
            Stats = dataset(SS,'ObsNames',obssnames','VarNames',{'SumofSquares'});
            Stats.DF = DF;
            Stats.MeanSquares = (SS./DF);
            Stats.Fvalue = F;           
            Stats.Pvalue = P;
            Stats.R = R;
            Stats.Direction = Direction;
            Stats.Intercept = Intercept;
        else
            New = dataset(SS(2,1),'ObsNames',char(VariableNames(vari)),'Varnames',{'SumofSquares'});
            New.DF = DF(2,1);
            New.MeanSquares = (SS(2,1)./DF(2,1));
            New.Fvalue = F(2,1);
            New.Pvalue = P(2,1);
            New.R = R(2,1);
            New.Direction = Direction(2,1);
            New.Intercept = Intercept(2,1);
            Stats(char(VariableNames(vari)),:) = New;
        end
        clear SS DF F P R  TargetVariable New 
        Mdls.(genvarname(char(EnsemblesNames(ensem)))).(genvarname(char(VariableNames(vari)))) = mdl2;
    end
    Results.(genvarname(char(EnsemblesNames(ensem)))) = Stats;
    Residuals.(genvarname(char(EnsemblesNames(ensem)))) = Yresid;
    
    clearvars -except Results EnsemblesNames VariableNames VariList run_max Residuals Mdls
end
save('ResultsWaterComparison','Results','Residuals','Mdls');
end

function [VariablesWater,EnsemblesWater] = CollateWaterData
Output_file = 'Results_water_All_Full.mat';
load(Output_file)
tst = Points.Models(:,2:10);
test = double(tst);
for i = 1: length(test)
    nrmodels(i,1) = length(find(isnan(test(i,:))~=1));
end
Points.Models.NrofModels = nrmodels;
EnsemblesValues = double(Points.Models(:,13:24));
EnsembleVariation = nanstd(EnsemblesValues,0,2)./sqrt(nrmodels);
Points.Models.SEMamongModels = EnsembleVariation;
Points.Comparison = Points.Models(:,[13:24,26,27]);
save(Output_file,'ResultsCombi','ResultsWeights','Parameters','VariationWeights','Points');
VariablesWater = dataset('File','Water extraction parameters.csv',...
   'Delimiter',',','ReadObsNames',true,'ReadVarNames',true);
EnsemblesWater = dataset(Points.Comparison,'ObsNames', Points.Models.Datapoint_name);
EnsemblesWater.CVamongEns = EnsemblesWater.SEMamongEnsembles.*sqrt(EnsemblesWater.NrofModels)./nanmean(double(EnsemblesWater(:,3:14)),2);
EnsemblesWater.CVamongModels = EnsemblesWater.STD_among_Models./ EnsemblesWater.Mean_among_Models;
save('TheWaterDataperKM','EnsemblesWater','VariablesWater')
end

function TarVar = normaVar(TarVar)
x_range_perc_low = prctile(TarVar(:,1),2.5);
x_range_norm = TarVar(:,1) - x_range_perc_low;
x_range_norm(x_range_norm<0) = 0;
x_range_perc = prctile(x_range_norm,97.5);
x_range = (x_range_norm./x_range_perc);
x_range(x_range>1) = 1;
TarVar = x_range;
clear x_range*
end