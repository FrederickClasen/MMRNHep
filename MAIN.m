%% FLUX SIMULATIONS PRESENTED IN CLASEN ET AL              %%%%%%%%%%%%%
% COBRA and RAVEN has to be added to the path to run the simulations

clear;
changeCobraSolver('mosek','all'); % SEE MOSEK WEBSITE FOR MORE DETAIL
setRavenSolver('mosek');          % SEE MOSEK WEBSITE FOR MORE DETAIL

%%%%%%%%%%%%%%          IMPORT GENERIC HEPATIC MODEL          %%%%%%%%%%%%%%%%%%%%%%

model = importExcelModel('data/models/xlsx/genericLiver.xlsx',false); % GENERIC HEPATIC MODEL MODEL

rxnsToDelete = {'AATAi','HMR_6404','PSERT'}; % BAD RXNS IN MODEL
model = removeReactions(model,rxnsToDelete,true,true,true); 

% RER AND OBJECTIVES
RERRxns = {'EXC_IN_C00007[s]','EXC_OUT_C00011[s]'};     
objective = {'MMRN_Biomass'};                           % OBJECTIVE FUNCTION
nonObjective = {'HMR_biomass_Renalcancer(with ATP)',... % RXNS THAT SHOULD BE BLOCKED AS POTENTIAL ADDITIONAL OBJECTIVES
                'EXC_OUT_ATP',...
                'EXC_OUT_Lipid_pool_biomass',...
                'EXC_OUT_Protein_pool_biomass',...
                'EXC_OUT_Nucleotide_pool_biomass',...
                'EXC_OUT_cofactors_vitamins',...
                'EXC_OUT_glycogen'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RUN FLUX SIMULATIONS FOR GENERIC AND CONSTRAINT-BASED GSMMs %% 

% load diet information
experiment = 'NORMAL';
load('data/diet.mat');
conditions = {'nonDEN_Liver_CD','nonDEN_Liver_WD','DEN_Liver_CD','DEN_AdjLiver_WD','DEN_Tumour_WD'};
fluxVector = zeros(length(model.rxns),(length(conditions)*2)+2); % FLUX VECTOR FOR TWO GENERIC MODELS AND CONDITION MODELS

counter = 0;
diets = {'WD','CD'};
% SIMULATE FOR TWO DIETS
for i=1:length(diets)
    counter = counter + 1;
    modelNames{counter} = char("generic"+string(diets{i}));
    
    dModel = model;   % generic model
    dModel = setParam(dModel,'ub',dModel.rxns(contains(dModel.rxns,'_IN_')),0);
    dModel = setParam(dModel,'ub',dModel.rxns(contains(dModel.rxns,'_BOTH_')),0);
    
    %SET DIET CONSTRAINT - THIS CAN BE CHANGED TO SWAP INDIVIDUAL
    %COMPONENTS WHEN RUNNING ONE DIET
    dModel = setParam(dModel,'ub',diet.aa.rxns,diet.aa.ub(:,i));
    dModel = setParam(dModel,'ub',diet.carbs.rxns,diet.carbs.ub(:,i));
    dModel = setParam(dModel,'ub',diet.lipids.rxns,diet.lipids.ub(:,i));
    dModel = setParam(dModel,'ub',diet.other.rxns,diet.other.ub(:,i));
    
    % SET RER AND OBJECTIVE
    dModel = setParam(dModel,'lb',RERRxns,[17,13]);
    dModel = setParam(dModel,'ub',RERRxns,[1000,1000]);
    dModel = setParam(dModel,'obj',objective,1);
    dModel = setParam(dModel,'eq',nonObjective,0);
    fba = solveLP(dModel,1);
    fluxVector(:,counter) = fba.x;
    
    for j=1:length(conditions)
        counter = counter + 1;
        modelNames{counter} = char(string(diets(i)) + '_' + string(conditions(j)));
        %model = models(i).model;
        fc = 'data/Eflux/' + string(conditions(j)) + '.csv';
        constraints = importdata(fc);
        csModel = dModel;
        csModel = addConstraints(csModel,constraints);  % EFLUX CONSTRAINT
        
        if contains(conditions(j),'nonDEN')
            csModel = setParam(csModel,'lb',RERRxns,[15,11]); % RER CONSTRAINT
            csModel = setParam(csModel,'ub',RERRxns,[1000,1000]);
        else
            csModel = setParam(csModel,'lb',RERRxns,[17,13]); % RER CONSTRAINT
            csModel = setParam(csModel,'ub',RERRxns,[1000,1000]);
        end
        
        csModel = setParam(csModel,'obj',objective,1);
        csModel = setParam(csModel,'ub',objective,1000);
        csModel = setParam(csModel,'eq',nonObjective,0);
        
        fba = solveLP(csModel,1);
        disp(fba.f);
        
        fluxVector(:,counter) = fba.x;
        
    end
    
end

tempModel = model;
T = table(tempModel.rxns,...
       constructEquations(tempModel,tempModel.rxns,true),...
       tempModel.grRules,...
       tempModel.subSystems);
T2 = array2table(fluxVector);
T = [T T2];
colheads = {'ID','EQUATION','GENEASSOCIATION','SUBSYSTEM'};
colheads = [colheads modelNames];
T.Properties.VariableNames = colheads;
writetable(T,'FLUXES.xlsx','Sheet',experiment,'WriteVariableNames',true);

