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
    {'SEMamongEnsembles'},{'STD_among_Models'},{'CVamongEns'},{'CVamongModels'}];
VariableNames = [{'RainfallAnnual'},{'PETAnnual'},{'PopSizeLog'},{'AgriculuturalPerc'},{'GrasslandPerc'},...
    {'PeatPerc'},{'UrbanPerc'},{'WoodlandPerc'},{'TemperatureMeanAnnual'},{'TemperatureMin'},...
    {'TemperatureSeasonality'},{'TemperatureRange'},{'RainfallSeasonality'},{'NrModels'},{'DEMMean'},{'SlopesMean'}];
VariList = [1:length(VariableNames)];
for ensem = 1:length(EnsemblesNames)
    load('TheWaterDataperKM.mat')
    TargetEnsemble = eval( ['EnsemblesWater.',char(EnsemblesNames(ensem))]);
    for vari = VariList
        str = sprintf('      Running variable %s in Ensemble %s',char(VariableNames(vari)),char(EnsemblesNames(ensem)));
        clc
        disp(str)
        obssnames = [{'Longitude'},{'Lattitude'},{'Nr_models'},{'LongxLat'},VariableNames(vari),{'Error'}];
        TargetVariable =  eval( ['VariablesWater.',char(VariableNames(vari))]);
        if min(TargetVariable)<0
            TargetVariable = TargetVariable+min(TargetVariable);
        end
        if min(TargetVariable)<0
            TargetVariable = TargetVariable+min(TargetVariable);
        end
        TarVar = normaVar(TargetVariable);
        LongitudeNorm = normaVar(VariablesWater.Longitude);
        LattitudeNorm = normaVar(VariablesWater.Lattitude);
        
        Y = TargetEnsemble;
        Longitude = LongitudeNorm;
        Lattitude = LattitudeNorm;
        NrofModels = VariablesWater.NrModels;
        
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
        XTarget = TargetVariable;
        ds2 = dataset(Yresid,XTarget);
        mdl2 = LinearModel.fit(ds2,'Yresid ~ XTarget');
        tbl2 = anova(mdl2,'component');
        
        SS(:,1)= [double(table(1,1));double(table(2,1));double(table(3,1));double(table(4,1));double(tbl2(1,1));double(tbl2(2,1))]; %#ok<*SAGROW>
        DF(:,1)=[double(table(1,2));double(table(2,2));double(table(3,2));double(table(4,2));double(tbl2(1,2));double(tbl2(2,2))];
        F(:,1)=[double(table(1,4));double(table(2,4));double(table(3,4));double(table(4,4));double(tbl2(1,4));double(tbl2(2,4))];
        P(:,1)=[double(table(1,5));double(table(2,5));double(table(3,5));double(table(4,5));double(tbl2(1,5));double(tbl2(2,5))];
        R(:,1) = [R_long;R_lat;R_Models;0;(min(max(mdl2.Rsquared.Adjusted,0),1));0];
        Direction(:,1)= [mdl.Coefficients.Estimate(2,1);mdl.Coefficients.Estimate(3,1);mdl.Coefficients.Estimate(4,1);...
            mdl.Coefficients.Estimate(5,1);mdl2.Coefficients.Estimate(2,1);0];
        %end
        if vari == VariList(1)
            Stats = dataset(SS,'ObsNames',obssnames','VarNames',{'SumofSquares'});
            Stats.DF = DF;
            Stats.MeanSquares = (SS./DF);
            Stats.Fvalue = F;
            Stats.Pvalue = P;
            Stats.R = R;
            Stats.Direction = Direction;
        else
            New = dataset(SS(5,1),'ObsNames',char(VariableNames(vari)),'Varnames',{'SumofSquares'});
            New.DF = DF(5,1);
            New.MeanSquares = (SS(5,1)./DF(5,1));
            New.Fvalue = F(5,1);
            New.Pvalue = P(5,1);
            New.R = R(5,1);
            New.Direction = Direction(5,1);
            Stats(char(VariableNames(vari)),:) = New;
        end
        clear SS DF F P R  TargetVariable New
    end
    Results.(genvarname(char(EnsemblesNames(ensem)))) = Stats;
    clearvars -except Results EnsemblesNames VariableNames VariList run_max
end
save('ResultsWaterComparison','Results');
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