
%%
addpath(genpath('/Users/clasenf/OneDrive - The Francis Crick Institute/PhD Pathway/Mouse model paper/revision/MMRNHep'));

%MMRN_SBML = readCbModel('data/models/sbml/MMRN.xml');

% 
% MMRN_EXCEL = setParam(MMRN_EXCEL,'obj','HMR_biomass_Renalcancer',1);
% obj = MMRN.rxns(find(MMRN.c == 1));
% solveLP(MMRN_EXCEL)

% for i=1:length(MMRN_EXCEL.mets)
%     disp(i)
%     model = addExchangeRxns(MMRN_EXCEL,'both',MMRN_EXCEL.mets(i));
%     fba = solveLP(model);
%     if abs(fba.f) > 0.1
%         disp(MMRN_EXCEL.mets(i));
%         break
%     end
% end
%[MMRN_EXCEL] = generateRules(MMRN_EXCEL);
%writeCbToSBML(MMRN_EXCEL,'data/models/sbml/MMRN.xml');

%MMRN = importExcelModel('data/models/xlsx/MMRN.xlsx');
%MMRN = setParam(MMRN,'obj','MMRN10785',1);
%solveLP(MMRN)
%model = extractSubNetwork(MMRN_EXCEL, MMRN_EXCEL.rxns(SConsistentRxnBool))
%[MMRN] = generateRules(MMRN);

MMRN = importExcelModel('data/models/xlsx/MMRN_revised.xlsx');
MMRN = setParam(MMRN,'obj','MMRNR10781',1);
solveLP(MMRN);
% for i=1:length(MMRN.mets)
%     disp(i)
%     model = addExchangeRxns(MMRN,'both',MMRN.mets(i));
%     fba = solveLP(model);
%     if abs(fba.f) > 0.1
%         disp(MMRN.mets(i));
%         break
%     end
% end

%% ADD ANNOTATIONS

clear
MMRN = importExcelModel('data/models/xlsx/MMRN_revised.xlsx');

opts = spreadsheetImportOptions("NumVariables", 24);

% Specify sheet and range
opts.Sheet = "RXNS";
opts.DataRange = "A2:X10832";

% Specify column names and types
opts.VariableNames = ["VarName1", "ID", "NAME", "EQUATION", "GENEASSOCIATION", "LOWERBOUND", "UPPERBOUND", "OBJECTIVE", "COMPARTMENT", "MIRIAM", "SUBSYSTEM", "REPLACEMENTID", "REFERENCE", "CONFIDENCESCORE", "rxnKEGGID", "rxnMetaNetXID", "rxnREACTOMEID", "rxnRheaID", "rxnRheaMasterID", "rxnBiGGID", "rxnRheaFinal", "rxnSBOTerms", "rxnECNumber", "RetiredID"];
opts.VariableTypes = ["char", "char", "char", "char", "categorical", "char", "char", "char", "char", "char", "categorical", "char", "char", "double", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "ID", "NAME", "EQUATION", "LOWERBOUND", "UPPERBOUND", "OBJECTIVE", "COMPARTMENT", "MIRIAM", "REPLACEMENTID", "REFERENCE", "rxnKEGGID", "rxnMetaNetXID", "rxnREACTOMEID", "rxnRheaID", "rxnRheaMasterID", "rxnBiGGID", "rxnRheaFinal", "rxnSBOTerms", "rxnECNumber", "RetiredID"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "ID", "NAME", "EQUATION", "GENEASSOCIATION", "LOWERBOUND", "UPPERBOUND", "OBJECTIVE", "COMPARTMENT", "MIRIAM", "SUBSYSTEM", "REPLACEMENTID", "REFERENCE", "rxnKEGGID", "rxnMetaNetXID", "rxnREACTOMEID", "rxnRheaID", "rxnRheaMasterID", "rxnBiGGID", "rxnRheaFinal", "rxnSBOTerms", "rxnECNumber", "RetiredID"], "EmptyFieldRule", "auto");

% Import the data
rxns = readtable("/Users/clasenf/OneDrive - The Francis Crick Institute/PhD Pathway/Mouse model paper/revision/MMRNHep/data/models/xlsx/MMRN_revised.xlsx", opts, "UseExcel", false);

MMRN.rxnKEGGID = cellstr(rxns.rxnKEGGID);
MMRN.rxnBiGGID = cellstr(rxns.rxnBiGGID);
MMRN.rxnMetaNetXID = cellstr(rxns.rxnMetaNetXID);
MMRN.rxnREACTOMEID = cellstr(rxns.rxnREACTOMEID);
MMRN.rxnRheaID = cellstr(rxns.rxnRheaFinal);
MMRN.rxnSBOTerms = cellstr(rxns.rxnSBOTerms);
MMRN.rxnECNumbers = cellstr(rxns.rxnECNumber);

opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = "METS";
opts.DataRange = "A2:P7188";

% Specify column names and types
opts.VariableNames = ["VarName1", "ID", "NAME", "UNCONSTRAINED", "MIRIAM", "COMPOSITION", "InChI", "COMPARTMENT", "REPLACEMENTID", "CHARGE", "hmdb", "kegg", "pubchem", "chebi", "lipidmaps", "metSBOTerms"];
opts.VariableTypes = ["char", "char", "categorical", "char", "categorical", "categorical", "char", "categorical", "char", "char", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "ID", "UNCONSTRAINED", "InChI", "REPLACEMENTID", "CHARGE", "hmdb", "kegg", "pubchem", "chebi", "lipidmaps", "metSBOTerms"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "ID", "NAME", "UNCONSTRAINED", "MIRIAM", "COMPOSITION", "InChI", "COMPARTMENT", "REPLACEMENTID", "CHARGE", "hmdb", "kegg", "pubchem", "chebi", "lipidmaps", "metSBOTerms"], "EmptyFieldRule", "auto");

% Import the data
mets = readtable("/Users/clasenf/OneDrive - The Francis Crick Institute/PhD Pathway/Mouse model paper/revision/MMRNHep/data/models/xlsx/MMRN_revised.xlsx", opts, "UseExcel", false);

MMRN.metKEGGID = cellstr(mets.kegg);
MMRN.metHMDBID = cellstr(mets.hmdb);
MMRN.metPubChemID = cellstr(mets.pubchem);
MMRN.metChEBIID = cellstr(mets.chebi);
MMRN.metLIPIDMAPSID = cellstr(mets.lipidmaps);
MMRN.metSBOTerms = cellstr(mets.metSBOTerms);

opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "GENES";
opts.DataRange = "A2:K3467";

% Specify column names and types
opts.VariableNames = ["VarName1", "NAME", "MIRIAM", "SHORTNAME", "COMPARTMENT", "ensembl_gene", "ensembl_transcript", "ensembl_protein", "entrez_gene", "refseq_gene", "uniprot_gene"];
opts.VariableTypes = ["char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["VarName1", "NAME", "MIRIAM", "SHORTNAME", "COMPARTMENT", "ensembl_gene", "ensembl_transcript", "ensembl_protein", "entrez_gene", "refseq_gene", "uniprot_gene"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "NAME", "MIRIAM", "SHORTNAME", "COMPARTMENT", "ensembl_gene", "ensembl_transcript", "ensembl_protein", "entrez_gene", "refseq_gene", "uniprot_gene"], "EmptyFieldRule", "auto");

% Import the data
genes = readtable("/Users/clasenf/OneDrive - The Francis Crick Institute/PhD Pathway/Mouse model paper/revision/MMRNHep/data/models/xlsx/MMRN_revised.xlsx", opts, "UseExcel", false);

MMRN.geneEnsemblGeneID = cellstr(genes.ensembl_gene);
MMRN.geneEnsemblTranscriptID = cellstr(genes.ensembl_transcript);
MMRN.geneEnsemblProteinID = cellstr(genes.ensembl_protein);
MMRN.geneEntrezID = cellstr(genes.entrez_gene);
MMRN.geneRefseqID = cellstr(genes.refseq_gene);
MMRN.geneUniprotID = cellstr(genes.uniprot_gene);

MMRN = generateRules(MMRN);

%%
model = MMRN;
[SConsistentMetBool, SConsistentRxnBool, SInConsistentMetBool, SInConsistentRxnBool, unknownSConsistencyMetBool, unknownSConsistencyRxnBool, model] = findStoichConsistentSubset(model);
MMRN_consistent = removeRxns(model,MMRN.rxns(SInConsistentRxnBool));
MMRN_consistent = removeRxns(model,MMRN.rxns(unknownSConsistencyRxnBool));
%writeCbToSBML(MMRN_consistent,'data/models/sbml/MMRN_consistent.xml');
%model.metNames(SInConsistentMetBool)

%%

MMRN = setParam(MMRN,'obj','MMRNR10760',1);
solveLP(MMRN)
MMRN_consistent = setParam(MMRN_consistent,'obj','MMRNR10760',1);
solveLP(MMRN_consistent)

%%
writeCbModel(MMRN,'data/models/sbml/MMRN.xml');
%writeCbToSBML(MMRN,'data/models/sbml/MMRN.xml');
%writeCbToSBML(MMRN_consistent,'data/models/sbml/MMRN_consistent.xml');

%%

exportToExcelFormat(MMRN_consistent,'data/models/xlsx/MMRN_consistent.xlsx');
%MMRN = setParam(MMRN,'obj','MMRNR10781',1);
%MMRN_consistent = setParam(MMRN_consistent,'obj','MMRNR10781',1);
%MMRN_fba = solveLP(MMRN,1)
%MMRN_consistent_fba = solveLP(MMRN_consistent,1)

%%

Human1 = readCbModel('Human-GEM.xml');
[SConsistentMetBool, SConsistentRxnBool, SInConsistentMetBool, SInConsistentRxnBool, unknownSConsistencyMetBool, unknownSConsistencyRxnBool, model] = findStoichConsistentSubset(Human1);
model.metNames(find(SInConsistentMetBool))

%%

cd = importExcelModel('data/models/xlsx/genericLiverLFD.xlsx',false);
cd = generateRules(cd);
solveLP(cd)
wd = importExcelModel('data/models/xlsx/genericLiverWD.xlsx',false);
wd = generateRules(wd);
solveLP(wd)
writeCbToSBML(cd,'data/models/sbml/MMRNHep-CD.xml');
writeCbToSBML(wd,'data/models/sbml/MMRNHep-WD.xml');





















