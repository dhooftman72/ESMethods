function Ensembles_Analyses_UK(validation_set,testset,runMax)
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
for validation_set = validation_set % To be added with function call
    clearvars -except validation_set testset runMax
    Parameters.testRun = testset;
    Parameters.runMax = runMax;
    Parameters.Ensemble_Names = {'Best Model Deviance';'Best Model Rho';... % Reference (1,2)
        'Mean';'Median';'PCA_weighting';'Median_among';'RegressAmong';'RegressAmongFlipped';... % Non-Informed (3-8)
        'MaxentDev';'MaxentRho';'MaxentDevFlip';'MaxentRhoFlip';'CorCoef';'CorrCoefFlip';'UniqueUpWeight';'UniqueDownWeight';... % Non-Informed (9-16)
        'HalfInformed_Deviance'; 'HalfInformed_Rho';'HalfInformed_BestDev';'HalfInformed_BestRho';...% Half-Informed (17 (cut-off) - 20)
        'HalfInformed_Bagging';'HalfInformed_IterBaggingDev'; 'HalfInformed_IterBaggingRho'}; % Half-Informed (22 - 23)
    
    % set all model, comparator values and further parameters
    [Parameters, Models,Comparator] = DefintionSet(validation_set,Parameters); 
    Parameters = ClusterTest(Parameters);
    if Parameters.testRun ~= 1 && Parameters.testRun ~= 2  
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
    [ResultsCombi,ResultsWeights] = JoinFuncUK(Parameters); %#ok<NASGU,ASGLU>
    Output_file = [Parameters.output_file,'_Combined'];
    save(Output_file,'ResultsCombi','ResultsWeights','Parameters');
    str = sprintf('Ready with validation set %s ',char(Parameters.SetNames(validation_set)));
    display ('  ')
    display ('  ')   
    disp(str)  
    display ('Combining all Runs')
end
display ('  ')
display ('  ')   
display ('Ready with full program')
end