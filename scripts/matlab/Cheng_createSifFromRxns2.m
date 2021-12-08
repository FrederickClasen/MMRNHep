function Cheng_createSifFromRxns(model,rxns,out,compartments,use_rev,fileName)
% Cheng_createSifFromRxns(model,rxns,out,compartments,use_rev,fileName)
% function createSifFromRxns
%
% This function creates a Cytoscape-readable .sif file for
% metabolite-gene/reaction networks. Only a GEM must be specified. If run
% on default, the resulting file is a non-compartmentalized metabolite-gene 
% non-directional network of all genes found in the model.
% The output is saved in the current directory with the name
% MetsToGenes.sif'.
%
% INPUT
% model:    RAVEN-format GEM
% rxns:     (optional) A cell array of strings containing the list of
%           reactions for which the network should be reconstructed. Default is all
%           reactions in the model. The cell array of strings must match the .rxns
%           field in the model to be included in the file.
% out:      (optional) A string, either 'rm','gm' (default) or 'rg' whether a
%           metabolite-reaction network, a metabolite-gene network, or a gene-reaction
%           network should be created.
% compartments: (optional) a logical, indicating whether metabolites in the
%               file should be considered unique if belonging to different compartments.
%               Default is false.
% use_rev:  (optional) A logical, indicating whether reversible reactions are 
%           written in the file in both directions. Default is false.
% fid:      (optional) A file ID for the output file. Default is
%           "MetsToGenes.sif"
function Cheng_createSifFromRxns(model,rxns,out,compartments,use_rev,fileName)
% Cheng_createSifFromRxns(model,rxns,out,compartments,use_rev,fileName)
% function createSifFromRxns
%
% This function creates a Cytoscape-readable .sif file for
% metabolite-gene/reaction networks. Only a GEM must be specified. If run
% on default, the resulting file is a non-compartmentalized metabolite-gene 
% non-directional network of all genes found in the model.
% The output is saved in the current directory with the name
% MetsToGenes.sif'.
%
% INPUT
% model:    RAVEN-format GEM
% rxns:     (optional) A cell array of strings containing the list of
%           reactions for which the network should be reconstructed. Default is all
%           reactions in the model. The cell array of strings must match the .rxns
%           field in the model to be included in the file.
% out:      (optional) A string, either 'rm','gm' (default) or 'rg' whether a
%           metabolite-reaction network, a metabolite-gene network, or a gene-reaction
%           network should be created.
% compartments: (optional) a logical, indicating whether metabolites in the
%               file should be considered unique if belonging to different compartments.
%               Default is false.
% use_rev:  (optional) A logical, indicating whether reversible reactions are 
%           written in the file in both directions. Default is false.
% fid:      (optional) A file ID for the output file. Default is
%           "MetsToGenes.sif"


model=importExcelModel('iHepatocytes2322_v4.xlsx',false);
rxns=table2array(cirrhosiscombined2(:,'RXNID'));

%u can use several outputs
Cheng_createSifFromRxns(model,rxns,'rm',true,true,'MetsToRXNS.sif')
Cheng_createSifFromRxns(model,rxns,'gm',true,false,'MetsToGenes.sif')
Cheng_createSifFromRxns(model,rxns,'rg',true,false,'RXNSToGenes.sif')

%for final network use
Cheng_createSifFromRxns(model,rxns,'rm',true,true,'steatosis.sif')