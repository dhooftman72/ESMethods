function   [Parameters, Models,Comparator] = DefintionSet(validation_set,Parameters)

%% Parameters various needed
Parameters.CutOffType = 17;
Parameters.max_its = 250;
Parameters.delta = [0.1;0.25];
Parameters.testRun = 0;
Parameters.runMax = 250;
Parameters.TProp = 0.5;
Parameters.Precision = 0.00001;
Parameters.NrCat = 4;
Parameters.ImprovementNrTake = 50;
Parameters.TimeOut = 120;

Parameters.Precision(2) = Parameters.Precision;%.*10;
Parameters.Precision = (1./Parameters.Precision);
Parameters.make_log = 0;
Parameters.ensemble = 0;
 
if Parameters.testRun == 1
    Parameters.runMax = 10;
    Parameters.TimeOut = 5;
    Parameters.max_its =25;
    Parameters.ImprovementNrTake = 5;
end
%% load data
if validation_set ==1
     Parameters.ServiceName = 'Carbon_tons';
    %GB set
    cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 0;
    Parameters.output_file = 'Results_Carbon_GB';
    Parameters.SetNames = {'Aries';'Invest_L2015';'transfer_L2015';'Henrys';'Baredo';'Kindermann';'Density';'NEInventory';'LPJveg';'Luci'};
        Models.Service = [Carbon.GB.(genvarname(char(Parameters.SetNames(1)))), Carbon.GB.(genvarname(char(Parameters.SetNames(2)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(3)))),Carbon.GB.(genvarname(char(Parameters.SetNames(4)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(5)))),Carbon.GB.(genvarname(char(Parameters.SetNames(6)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(7)))),Carbon.GB.(genvarname(char(Parameters.SetNames(8)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(9)))),Carbon.GB.(genvarname(char(Parameters.SetNames(10))))];
    Models.Areas = [Sizes.GB.(genvarname(char(Parameters.SetNames(1)))), Sizes.GB.(genvarname(char(Parameters.SetNames(2)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(3)))),Sizes.GB.(genvarname(char(Parameters.SetNames(4)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(5)))),Sizes.GB.(genvarname(char(Parameters.SetNames(6)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(7)))),Sizes.GB.(genvarname(char(Parameters.SetNames(8)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(9)))),Sizes.GB.(genvarname(char(Parameters.SetNames(10))))];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Carbon.GB.Invest_L2007),Parameters.data_set_max);
    Parameters.names = Carbon.GB.Forest;
    Comparator.Service = Carbon.GB.ForestResearch;
    Parameters.Sizes = Sizes.GB.ForestResearchArea;
    Parameters.UniqueMatrix =  Carbon.UniquenessMatrix;
elseif validation_set ==2
    Parameters.ServiceName = 'Water_flows';
    %All set
     cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 1;
    Parameters.output_file  = 'Results_water_All';
    Parameters.SetNames= {'Invest_L2015';'WaterWorld';'Growth_days';'transfer_L2015';'G2G';....
        'Aquaduct';'Decipher';'LPJ';'Luci'};
    Models.Service = [Flows.All.(genvarname(char(Parameters.SetNames(1)))), Flows.All.(genvarname(char(Parameters.SetNames(2)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(3)))),((Flows.All.(genvarname(char(Parameters.SetNames(4)))))./Sizes.All.Population),...
        Flows.All.(genvarname(char(Parameters.SetNames(5)))),Flows.All.(genvarname(char(Parameters.SetNames(6)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(7)))),Flows.All.(genvarname(char(Parameters.SetNames(8)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(9))))];
    Models.Areas = [Sizes.All.(genvarname(char(Parameters.SetNames(1)))), Sizes.All.(genvarname(char(Parameters.SetNames(2)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(3)))),Sizes.All.(genvarname(char(Parameters.SetNames(4)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(5)))),Sizes.All.(genvarname(char(Parameters.SetNames(6)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(7)))),Sizes.All.(genvarname(char(Parameters.SetNames(8)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(9))))];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Flows.All.Invest_L2015),data_set_max);
    Parameters.names = Flows.All.Catchment;
    Comparator.Service = Flows.All.TwentyYearFlow;
    Parameters.Sizes = Sizes.All.NFRASize;
    Parameters.UniqueMatrix =  Flows.UniquenessMatrix;
elseif validation_set == 3
     Parameters.ServiceName = 'Carbon_tons';
    %GB set
    cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 0;
    Parameters.output_file = 'Results_Carbon_GB_Meta';
    Parameters.SetNames = {'Aries';'InvestMeta';'TransferMeta';'Henrys';'Baredo';'Kindermann';'Density';'NEInventory';'LPJveg';'Luci'};
        Models.Service = [Carbon.GB.(genvarname(char(Parameters.SetNames(1)))), Carbon.GB.(genvarname(char(Parameters.SetNames(2)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(3)))),Carbon.GB.(genvarname(char(Parameters.SetNames(4)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(5)))),Carbon.GB.(genvarname(char(Parameters.SetNames(6)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(7)))),Carbon.GB.(genvarname(char(Parameters.SetNames(8)))),...
        Carbon.GB.(genvarname(char(Parameters.SetNames(9)))),Carbon.GB.(genvarname(char(Parameters.SetNames(10))))];
    Models.Areas = [Sizes.GB.(genvarname(char(Parameters.SetNames(1)))), Sizes.GB.(genvarname(char(Parameters.SetNames(2)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(3)))),Sizes.GB.(genvarname(char(Parameters.SetNames(4)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(5)))),Sizes.GB.(genvarname(char(Parameters.SetNames(6)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(7)))),Sizes.GB.(genvarname(char(Parameters.SetNames(8)))),...
        Sizes.GB.(genvarname(char(Parameters.SetNames(9)))),Sizes.GB.(genvarname(char(Parameters.SetNames(10))))];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Carbon.GB.Invest_L2007),Parameters.data_set_max);
    Parameters.names = Carbon.GB.Forest;
    Comparator.Service = Carbon.GB.ForestResearch;
    Parameters.Sizes = Sizes.GB.ForestResearchArea;
    Parameters.UniqueMatrix =  Carbon.UniquenessMatrix;
elseif validation_set == 4
     Parameters.ServiceName = 'Water_flows';
    %All set
     cd('TheData')
    load(Parameters.ServiceName)
    cd ..
    Parameters.make_log = 1;
    Parameters.output_file  = 'Results_water_All_Meta';
    Parameters.SetNames= {'InvestMeta';'WaterWorld';'Growth_days';'TransferMeta';'G2G';....
        'Aquaduct';'Decipher';'LPJ';'Luci'};
    Models.Service = [Flows.All.(genvarname(char(Parameters.SetNames(1)))), Flows.All.(genvarname(char(Parameters.SetNames(2)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(3)))),((Flows.All.(genvarname(char(Parameters.SetNames(4)))))./Sizes.All.Population),...
        Flows.All.(genvarname(char(Parameters.SetNames(5)))),Flows.All.(genvarname(char(Parameters.SetNames(6)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(7)))),Flows.All.(genvarname(char(Parameters.SetNames(8)))),...
        Flows.All.(genvarname(char(Parameters.SetNames(9))))];
    Models.Areas = [Sizes.All.(genvarname(char(Parameters.SetNames(1)))), Sizes.All.(genvarname(char(Parameters.SetNames(2)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(3)))),Sizes.All.(genvarname(char(Parameters.SetNames(4)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(5)))),Sizes.All.(genvarname(char(Parameters.SetNames(6)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(7)))),Sizes.All.(genvarname(char(Parameters.SetNames(8)))),...
        Sizes.All.(genvarname(char(Parameters.SetNames(9))))];
    Parameters.data_set_max =  size(Models.Service,2);
    %Models.Areas = ones(length(Flows.All.Invest_L2015),data_set_max);
    Parameters.names = Flows.All.Catchment;
    Comparator.Service = Flows.All.TwentyYearFlow;
    Parameters.Sizes = Sizes.All.NFRASize;
    Parameters.UniqueMatrix =  Flows.UniquenessMatrix;
end
end