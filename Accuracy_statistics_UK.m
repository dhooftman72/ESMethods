% Statistics function for Hooftman et al. EnsemblES
function [Outputs] = Accuracy_statistics_UK(testArray,Parameters,Szs)
% clean data set
testArray(isinf(testArray)==1) = NaN;
orginal_order = 1:1:length(testArray(:,1));
orginal_order = orginal_order';
a1=find((isnan(testArray(:,1))==1));
b=find((isnan(testArray(:,2))==1));
c = find(Szs <2);% for carbon to remove small framents
all = [a1;b;c];
a(:,1) = unique(all);
testArray(a,:) = [];  %#ok<*FNDSB>
orginal_order(a,:) = [];
missing_order = a;
clear a
clear a1
clear b
clear c
clear all
Outputs.datapoints = size(testArray,1);
%% Make log if applicable, but not for Ensembles
clear testVar
testVar = testArray;
if  Parameters.make_log == 1
    testVar = log10(testArray+1);
end
%% Correlation stats: Spearman: Overall Accuracy
isnnx = length(find(isnan(testVar(:,1)))==1);
isnny = length(find(isnan(testVar(:,2)))==1);
if (exist('testVar','var') == 1) && (isempty(testVar)~= 1) &&  (isnnx ~= length(testVar(:,1))) &&  (isnny ~= length(testVar(:,2)))
    [Outputs.RHO,Outputs.PVAL] = corr(testVar(:,1),testVar(:,2),'type','Spearman');
    Outputs.PVAL = single(Outputs.PVAL);
    Outputs.RHO = (round(Outputs.RHO.*Parameters.Precision(2)))./Parameters.Precision(2);
    %% Inverse deviance against a 1:1 line
    clear x_range y_range
    % Mean double normalised deviation: !Winzorisation!
    % Two sides and on range!!
    if Parameters.ensemble == 0
        x_range_perc_low = prctile(testVar(:,1),2.5);
        y_range_perc_low = prctile(testVar(:,2),2.5);
        x_range_norm = testVar(:,1) - x_range_perc_low;
        y_range_norm = testVar(:,2) - y_range_perc_low;
        x_range_norm(x_range_norm<0) = 0;
        y_range_norm(y_range_norm<0) = 0;
        x_range_perc = prctile(x_range_norm,97.5);
        y_range_perc = prctile(y_range_norm,97.5);
        x_range = (x_range_norm./x_range_perc);
        y_range = (y_range_norm./y_range_perc);
        x_range(x_range>1) = 1;
        y_range(y_range>1) = 1;
    elseif Parameters.ensemble == 1
        x_range= testVar(:,1);
        y_range = testVar(:,2);
        x_range(x_range<0) = 0;
        y_range(y_range<0) = 0;
        x_range(x_range>1) = 1;
        y_range(y_range>1) = 1;
    end
    clear x_range_* y_range_*
    clear deviation_point
    Outputs.deviation_point= abs(y_range-x_range); %% Accuracy per point
    Outputs.mean_double_deviation = 1- ((sum(abs(y_range-x_range)))/Outputs.datapoints); %%Accuracy overall
    Outputs.LeastSquares = -sum((abs(y_range-x_range)).^2);
    Outputs.std_double_deviation = nanstd(1-(abs(y_range-x_range))); % Output for sensitivities
    Outputs.mean_double_deviation = (round( Outputs.mean_double_deviation.*Parameters.Precision(2)))./Parameters.Precision(2);
    Outputs.std_double_deviation  = (round(Outputs.std_double_deviation.*Parameters.Precision(2)))./Parameters.Precision(2);
    % %% recreate x and y ranges in orginal order, for individual model point
    % storage to create ensembles
    clear xes yes
Outputs.xes(orginal_order,1) = x_range;
Outputs.xes(missing_order,1) = NaN;
Outputs.yes(orginal_order,1) = y_range;
Outputs.yes(missing_order,1) = NaN;
Outputs.deviation_point(orginal_order,1) = Outputs.deviation_point;
Outputs.deviation_point(missing_order,1) = NaN;
else 
   Outputs.RHO = -9999;
end
end