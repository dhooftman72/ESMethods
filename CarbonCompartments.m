clc
clear all
warning off %#ok<*WNOFF>
load('FR_data.mat')
for country = 1:1:2
    if country == 1
        Compartments = CompartmentListEngland;
    elseif country == 2
        Compartments = CompartmentListScotland;
    end
    TotalLength = length(Compartments.FID);
    Carbon_list = dataset(-999,'Varnames',char('FID'));
    Wierd_list =  dataset(-999,'Varnames',char('FID'));
    wierd_count = 0;
    for x = 1:1:TotalLength % all subcompartments
        if x/10 == ceil(x/10)
            clc
            display(x)
        end
        carbon = 0;
        Carbon_list.FID(x,1) = Compartments.FID(x);
        Specs = [Compartments.Pri_Spec_code(x),...
            Compartments.Sec_Spec_code(x),...
            Compartments.Third_Spec_code(x)];
        Ages = [Compartments.PriAge(x),...
            Compartments.SecAge(x),...
            Compartments.ThirdAge(x)];
        Covers = [Compartments.PriCover(x),...
            Compartments.SecCover(x),...
            Compartments.ThirdCover(x)];
        Yields = [Compartments.PriYieldCode(x),...
            Compartments.SecYieldCode(x),...
            Compartments.ThirdYieldCode(x)];
        Area = Compartments.Hectare(x);
        for i = 1:1:3 % for three species layers
            Spec = Specs(i);
            Age = Ages(i);
            Cover = Covers(i);
            Yield = Yields(i);
            %Spec changes, indentified in earlier runs
            if strcmp(Spec,'JUN') == 1 % Juniper
                wierd_count = wierd_count + 1;
                Wierd_list.FID(wierd_count,1) = Compartments.FID(x);
                Wierd_list.layer(wierd_count,1) = i;
                Wierd_list.Spec(wierd_count,1) = Spec;
                Wierd_list.Age(wierd_count,1) = Age;
                 Wierd_list.SubCompCode(wierd_count,1) = Compartments.SubCompCode(x);
                Wierd_list.Solved(wierd_count,1) = {'solved'};
                Spec = {'XC'};
            end
            if strcmp(Spec,'XBN') == 1 % Non Native other broadleaves
                wierd_count = wierd_count + 1;
                Wierd_list.FID(wierd_count,1) = Compartments.FID(x);
                Wierd_list.layer(wierd_count,1) = i;
                Wierd_list.Spec(wierd_count,1) = Spec;
                Wierd_list.Age(wierd_count,1) = Age;
                  Wierd_list.SubCompCode(wierd_count,1) = Compartments.SubCompCode(x);
                Wierd_list.Solved(wierd_count,1) = {'solved'};
                Spec = {'XB'};
            end
            if  strcmp(Spec,'None') ~= 0
                carbon(i) = 0; %#ok<*SAGROW>
            else % start looking for the carbon value per hectare
                Lookup = find(strcmp(Spec, Full_Species_list.Code(:)));
                if isempty(Lookup) == 1 % the species is unknown
                    wierd_count = wierd_count + 1;
                    Wierd_list.FID(wierd_count,1) = Compartments.FID(x);
                    Wierd_list.layer(wierd_count,1) = i;
                    Wierd_list.Spec(wierd_count,1) = Spec;
                    Wierd_list.Age(wierd_count,1) = Age;
                     Wierd_list.SubCompCode(wierd_count,1) = Compartments.SubCompCode(x);
                    Wierd_list.Solved(wierd_count,1) = {'Unsolved, will be removed'};
                    carbon(i) = 0;
                else % Species is known
                    Model = Full_Species_list.Model(Lookup);
                    Listtmp =  find(strcmp(Model, Yield_list.Species(:))); % Right species
                    Yield_listtmp = Yield_list(Listtmp,:);  %#ok<*FNDSB>
                    Listtmp = find(Yield_listtmp.ManagmentCode==0); % No thinning
                    Yield_listtmp = Yield_listtmp(Listtmp,:);
                    % Correct yield code, note that with no yield code a
                    % broader array of values are used, becaus ethis step
                    Listtmp = find(Yield_listtmp.YieldCode==Yield); % Correct yield code
                    Yield_listtmp1 = Yield_listtmp(Listtmp,:);
                    clear Lookup
                    if isempty(Yield_listtmp1) ~= 1
                        Yield_listtmp = Yield_listtmp1;
                    end
                    clear Yield_listtmp1
                    if Age == 2019
                        carbon(i) = 0;
                    else
                        Age_cat = ceil(Age/5)*5;
                        if Age_cat > 1000
                            wierd_count = wierd_count + 1;
                            Wierd_list.FID(wierd_count,1) = Compartments.FID(x);
                            Wierd_list.layer(wierd_count,1) = i;
                            Wierd_list.Spec(wierd_count,1) = Spec;
                            Wierd_list.Age(wierd_count,1) = Age;
                            Wierd_list.SubCompCode(wierd_count,1) = Compartments.SubCompCode(x);
                            Wierd_list.Solved(wierd_count,1) = {'Unsolved because of Age, will be removed'};
                            Carbon_ha = 0;
                        else
                            if Age_cat > 200 % Correct for very old ages
                                Age_cat = 200;
                            end
                            if Age ~= 0
                                Listtmp = find(Yield_listtmp.EndYear == Age_cat); % Years
                                Yield_listtmp = Yield_listtmp(Listtmp,:);
                                Carbon_ha = nanmean(Yield_listtmp.CumCarbon);
                            else
                                Listtmp = find(Yield_listtmp.EndYear == 5);
                                Yield_listtmp = Yield_listtmp(Listtmp,:);
                                Carbon_ha = (nanmean(Yield_listtmp.CumCarbon))./5;
                            end
                        end
                        carbon(i)  = Carbon_ha.*(Cover/100).*Area;
                        if isnan(carbon(i))== 1
                            wierd_count = wierd_count + 1;
                            Wierd_list.FID(wierd_count,1) = Compartments.FID(x);
                            Wierd_list.layer(wierd_count,1) = i;
                            Wierd_list.Spec(wierd_count,1) = Spec;
                            Wierd_list.Age(wierd_count,1) = Age;
                              Wierd_list.SubCompCode(wierd_count,1) = Compartments.SubCompCode(x);
                            Wierd_list.Solved(wierd_count,1) = {'Unsolved unknown, will be removed'};
                            carbon(i) = 0;
                            
                        end
                    end % end age test
                    clear Yield_listtmp
                    clear Listtmp
                    clear Age_cat
                end % end known species
            end % species is known
            clear Spec
            clear Age
            clear Cover
            clear Yield
            clear Model
            clear i
        end % the three layers
        Carbon_list.CarbonLY1(x,1) = carbon(1);
        Carbon_list.CarbonLY2(x,1) = carbon(2);
        Carbon_list.CarbonLY3(x,1) = carbon(3);
        Carbon_list.TotalCarbon(x,1) = sum(carbon);
        Carbon_list.CarbonHa(x,1) = Carbon_list.TotalCarbon(x,1)./Area;
        Carbon_list.AreaHa(x,1) = Area;
        clear carbon
        clear Specs
        clear Ages
        clear Covers
        clear Yields
        clear Area
        clear Carbon_ha
        
    end % end of compartments
    clear x
    clear TotalLength
    clear wierd_count
    if country == 1
        Carbon_list_England = Carbon_list;
        Wierd_list_England = Wierd_list;
    elseif  country == 2
        Carbon_list_Scotland = Carbon_list;
        Wierd_list_Scotland = Wierd_list;
    end
clear Carbon_list
clear Wierd_list
end
clear country