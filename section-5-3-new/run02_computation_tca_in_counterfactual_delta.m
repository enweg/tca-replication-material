clear; clc;
addpath("./functions/")
addpath("/Applications/Dynare/6.3-arm64/matlab/")  % change this to your Dynare path
addpath("../tca-matlab-toolbox-v0.3.0")
addpath("../tca-matlab-toolbox-v0.3.0/models")

%% --- GLOBAL OPTIONS

runOptions = struct();
runOptions.maxHorizon = 12; 
runOptions.saveFileName = 'resultsDelta-from0p001-by0p001-to0p05.mat';

% --- Set the possible values for delta
deltaValues = 0.001:0.001:0.05;

%% --- BASELINE MODEL

% warmup run
cd dynare-models; dynare model_extensive_gap.mod; cd ..; addpath('dynare-models'); clc;
addpath('dynare-models');

inflationDecompositions = zeros(runOptions.maxHorizon + 1, 5, numel(deltaValues));

for i = 1:numel(deltaValues)
    global M_ options_ oo_;
    delta = deltaValues(i);
    set_param_value('delta', delta);
    options_.order = 1; 
    options_.irf = runOptions.maxHorizon; 
    options_.nograph = 1; 
    options_.noprint = 1;
    if isfield(oo_, 'irfs'); oo_.irfs = struct(); end
    steady; 
    [info, oo_] = stoch_simul(M_, options_, oo_, []);

    resultsBaseline = struct();
    resultsBaseline.M_ = M_; 
    resultsBaseline.options_ = options_; 
    resultsBaseline.oo_ = oo_; 

    resultsBaseline = computeEffects(runOptions, resultsBaseline, false); 

    irfs = resultsBaseline.irfsVarma; 
    effectConsumptionChannel = resultsBaseline.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue.effectsConsumption; 
    effectInvestmentChannel = resultsBaseline.TCAThroughConsumptionInvestmentInvariant.orderCBeforeITrue.effectsInvestment; 

    totalCons = irfs(1, 3, :); 
    totalCons = totalCons(:);
    totalInv = irfs(2, 3, :); 
    totalInv = totalInv(:);
    totalInfl = irfs(3, 3, :);
    totalInfl = totalInfl(:);
    consChannel = effectConsumptionChannel(3, 1, :); 
    consChannel = consChannel(:);
    invChannel = effectInvestmentChannel(3, 1, :);
    invChannel = invChannel(:); 

    inflationDecompositions(:, :, i) = [totalCons totalInv totalInfl consChannel invChannel];
end

results = struct();
results.inflationDecompostions = inflationDecompositions; 
results.columnNames = {'totalCons', 'totalInv', 'totalInfl', 'consChannel', 'invChannel'};
results.deltaValues = deltaValues;

save(fullfile('outputs/tca-in-counterfactual', runOptions.saveFileName), 'results');

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

