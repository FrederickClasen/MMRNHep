%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%          FLUX SIMULATIONS PRESENTED IN CLASEN ET AL              %%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%          ADD THE PATHS TO RAVEN AND COBRA          %%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('/Users/clasenf/OneDrive - The Francis Crick Institute/PhD Pathway/Mouse model paper/DataRepo/git')); % CURRENT PATH
addpath(genpath('/Users/clasenf/GEM/Raven/')); % RAVEN TOOLBOX
addpath(genpath('/Users/clasenf/GEM/COBRA/')); % COBRA TOOLBOX
changeCobraSolver('mosek','all');              % SEE MOSEK WEBSITE FOR MORE DETAIL
setRavenSolver('mosek');                       % SEE MOSEK WEBSITE FOR MORE DETAIL


%%%%%%%%%%%%%%          IMPORT GENERIC DIET CONSTRAINT GSMMS          %%%%%%%%%%%%%%%%%%%%%%

models(1).name = 'WD';
models(2).name = 'CD';
models(1).model = importExcelModel('data/models/xlsx/genericLiverWD.xlsx',false); % GENERIC WD MODEL
models(2).model = importExcelModel('data/models/xlsx/genericLiverCD.xlsx',false); % GENERIC CD MODEL

% RER AND OBJECTIVES
RERRxns = {'EXC_IN_C00007[s]','EXC_OUT_C00011[s]'};
objective = {'MMRN_Biomass'}; % OBJECTIVE FUNCTION
nonObjective = {'HMR_biomass_Renalcancer(with ATP)',... % RXNS THAT SHOULD BE BLOCKED AS POTENTIAL ADDITIONAL OBJECTIVES
                'EXC_OUT_ATP',...
                'EXC_OUT_Lipid_pool_biomass',...
                'EXC_OUT_Protein_pool_biomass',...
                'EXC_OUT_Nucleotide_pool_biomass',...
                'EXC_OUT_cofactors_vitamins',...
                'EXC_OUT_glycogen'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RUN FLUX SIMULATIONS FOR GENERIC AND CONSTRAINT-BASED GSMMs %%

conditions = {'nonDEN_Liver_CD','nonDEN_Liver_WD','DEN_Liver_CD','DEN_AdjLiver_WD','DEN_Tumour_WD'};
fluxVector = zeros(length(models(1).model.rxns),(length(conditions)*2)+2); % FLUX VECTOR FOR TWO GENERIC MODELS AND CONDITION MODELS

counter = 0;

for i=1:length(models)
    counter = counter + 1;
    modelNames{counter} = models(i).name;
    model = models(i).model;
    model = setParam(model,'lb',RERRxns,[17,13]);
    model = setParam(model,'obj',objective,1);
    model = setParam(model,'eq',nonObjective,0);
    fba = solveLP(model,1);
    fluxVector(:,counter) = fba.x;
    for j=1:length(conditions)
        counter = counter + 1;
        modelNames{counter} = char(string(models(i).name) + '_' + string(conditions(j)));
        model = models(i).model;
        fc = 'data/Eflux/' + string(conditions(j)) + '.csv';
        constraints = importdata(fc);
        outModel = addConstraints(model,constraints);  % EFLUX CONSTRAINT
        
        if contains(conditions(j),'nonDEN')
            outModel = setParam(outModel,'lb',RERRxns,[15,11]); % RER CONSTRAINT
        else
            outModel = setParam(outModel,'lb',RERRxns,[17,13]); % RER CONSTRAINT
        end
        
        outModel = setParam(outModel,'obj',objective,1);
        outModel = setParam(outModel,'ub',objective,1000);
        outModel = setParam(outModel,'eq',nonObjective,0);
        
        fba = solveLP(outModel,1);
        disp(fba.f);
        
        fluxVector(:,counter) = fba.x;
        
    end
    
end

tempModel = models(1).model;
T = table(tempModel.rxns,...
       constructEquations(tempModel,tempModel.rxns,true),...
       tempModel.grRules,...
       tempModel.subSystems);
T2 = array2table(fluxVector);
T = [T T2];
colheads = {'ID','EQUATION','GENEASSOCIATION','SUBSYSTEM'};
colheads = [colheads modelNames];
T.Properties.VariableNames = colheads;
writetable(T,'TableS3.xlsx','Sheet','TableS3','WriteVariableNames',true);











