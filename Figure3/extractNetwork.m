%% DEFINE PATHS AND INITIATE COBRA AND RAVEN %%



addpath(genpath('/Users/clasenf/OneDrive - The Francis Crick Institute/PhD Pathway/Mouse model paper/DataRepo/MMRNHep')); 
addpath(genpath('/Users/clasenf/GEM/Raven/'));                                     % RAVEN toolbox
addpath(genpath('/Users/clasenf/GEM/COBRA/'));                                     % COBRA toolbox
changeCobraSolver('mosek','all');
setRavenSolver('mosek');



%% EXTRACT SUBNETWORKS

%model = importExcelModel('../data/models/xlsx/genericLiverWD.xlsx',false);

rxns = importdata('Figure3DRxns.txt');
subModel = extractSubNetwork(model, rxns);
subModel=simplifyModel(subModel);
subModel=removeMets(subModel,{'H2O','CO2','O2','H+','HCO3-','Na+',...
                     'CoA','Pi','PPi','AMP','ADP','ATP',...
                     'CMP','CDP','CTP','GMP','GDP','GTP',...
                     'UTP','UDP','IDP','IMP','K+','NAD+',...
                     'NADH','NADP+','NADPH','PAP','PAPS',...
                     'FAD','FADH2','H2O2','O2-'},true);
Cheng_createSifFromRxns(subModel,rxns,'rm',true,true,'network/rxn2Mets.sif');
Cheng_createSifFromRxns(subModel,rxns,'rg',true,true,'network/rxn2Gene.sif');



