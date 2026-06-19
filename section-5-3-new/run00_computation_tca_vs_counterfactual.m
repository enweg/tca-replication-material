clear; clc;
addpath("./functions/")
addpath("/Applications/Dynare/6.3-arm64/matlab/")  % change this to your Dynare path
addpath("../tca-matlab-toolbox-v0.3.0")
addpath("../tca-matlab-toolbox-v0.3.0/models")

%% --- GLOBAL OPTIONS

runOptions = struct();
runOptions.maxHorizon = 12; 

%% --- BASELINE MODEL

cd dynare-models; dynare model_extensive_gap.mod; cd ..; clc; 
resultsBaseline = struct();
resultsBaseline.M_ = M_; 
resultsBaseline.options_ = options_; 
resultsBaseline.oo_ = oo_; 

resultsBaseline = computeEffects(runOptions, resultsBaseline, false); 

save('outputs/tca-vs-counterfactual/resultsBaseline.mat', 'resultsBaseline');

%% --- COUNTERFACTUAL SIGMA -> INFINITY

cd dynare-models; dynare model_extensive_gap_no_consumption.mod; cd ..; clc; 
resultsCounterfactualNoConsumption = struct();
resultsCounterfactualNoConsumption.M_ = M_; 
resultsCounterfactualNoConsumption.options_ = options_; 
resultsCounterfactualNoConsumption.oo_ = oo_; 

resultsCounterfactualNoConsumption = computeEffects(runOptions, resultsCounterfactualNoConsumption, true); 

save('outputs/tca-vs-counterfactual/resultsCounterfactualNoConsumption.mat', 'resultsCounterfactualNoConsumption');

%% --- COUNTERFACTUAL PHI_I -> INFINITY

cd dynare-models; dynare model_extensive_gap_no_investment.mod; cd ..; clc; 
resultsCounterfactualNoInvestment = struct();
resultsCounterfactualNoInvestment.M_ = M_; 
resultsCounterfactualNoInvestment.options_ = options_; 
resultsCounterfactualNoInvestment.oo_ = oo_; 

resultsCounterfactualNoInvestment = computeEffects(runOptions, resultsCounterfactualNoInvestment, true); 

save('outputs/tca-vs-counterfactual/resultsCounterfactualNoInvestment.mat', 'resultsCounterfactualNoInvestment');

%% --- FUNCTIONS

function results = computeEffects(runOptions, results, isCounterfactual)
    results = getTotalEffects(runOptions, results); 

    if ~isCounterfactual

        % order = (c, I, pi}
        runOptions.orderCBeforeI = true; 
        results = getTCAThroughConsumptionInvestmentSensitive(runOptions, results);
        results = getTCAThroughConsumptionInvestmentInvariant(runOptions, results); 
        results = getTCAClosestCounterfactual(runOptions, results); 

        % order = (I, c, pi)
        runOptions.orderCBeforeI = false; 
        results = getTCAThroughConsumptionInvestmentSensitive(runOptions, results);
        results = getTCAThroughConsumptionInvestmentInvariant(runOptions, results); 
        results = getTCAClosestCounterfactual(runOptions, results); 
    end
end

