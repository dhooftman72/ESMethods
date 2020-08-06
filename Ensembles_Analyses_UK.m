function Ensembles_Analyses_UK(validation_sets,testset,runMax)
% Clean and wipe
clc
close all hidden
if matlabpool('size') ~= 0
    matlabpool close
end
warning off
if exist('Output_Dir','dir') ~= 7
mkdir('Output_Dir')
end
%% join validation data with comparator
for validation_set = validation_sets % To be added with function call
    clearvars -except validation_set testset runMax
    Parameters.testRun = testset;
    Parameters.runMax = runMax;
    Parameters.Ensemble_Names = {'Best Model Deviance';'Best Model Rho';... % Reference (1,2)
        'Mean';'Median';'PCA_weighting';'Median_among';'RegressAmong';... % Non-Informed (3-7)
        'MaxentDev';'MaxentRho';'MaxentLeast';'CorCoef';'GridSize';'UniqueUpWeight';...% Non-Informed (8-13)
        'UniqueDownWeight';'AllInformed';... % Non-Informed (14-15)
        'HalfInformed_Deviance'; 'HalfInformed_Rho';'HalfInformed_BestDev';'HalfInformed_BestRho';...% Half-Informed (16 (cut-off) - 19)
        'HalfInformed_Bagging';'HalfInformed_IterBaggingDev'; 'HalfInformed_IterBaggingRho';'HalfInformed_IterBaggingLeast'}; % Half-Informed (20 - 23)
        %'RegressAmongFlipped';'MaxentDevFlip';'MaxentRhoFlip';
    
        % set all model, comparator values and further parameters
    [Parameters, Models,Comparator] = DefintionSet(validation_set,Parameters); 
    Parameters = ClusterTest(Parameters);
    if Parameters.testRun == 0
        display('Paralel job')
        job = createJob('configuration', 'Full');
        for run= 1:Parameters.runMax
            createTask(job, @TheRuns_UK, 0,{Parameters,Models,Comparator,run});
        end
        submit(job);
        waitForState(job, 'finished');
        %results = getAllOutputArguments(job);
        destroy(job)
    else
        display('No paralel')
        for run= 1:Parameters.runMax
           TheRuns_UK(Parameters,Models,Comparator,run);
        end
    end
    display ('  ')
    display ('  ')
    display ('Combining all Runs')
    [ResultsCombi,ResultsWeights,VariationWeights] = JoinFuncUK(Parameters); %#ok<NASGU,ASGLU>
    
    if Parameters.testRun ~= 2
        Output_file = [Parameters.output_file,'_Bootstrap'];
        save(Output_file,'ResultsCombi','ResultsWeights','Parameters','VariationWeights');
    else
         Output_file = [Parameters.output_file,'_Full'];
         cd('Output_Dir')
         load(Output_file);
         cd ..
         clear Results Weighting
         save(Output_file,'ResultsCombi','ResultsWeights','Parameters','VariationWeights','Points');
    end
    str = sprintf('Ready with validation set %s ',Parameters.output_file);
    display ('  ')
    display ('  ')   
    disp(str)  
end
display ('  ')
display ('  ')   
display ('Ready with full program')
end
