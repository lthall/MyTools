function [ cost ] = run_and_cost_function( variables )
%DIPOLE Summary of this function goes here
%   Detailed explanation goes here

global Opt_Data;
HFSS_Test = Opt_Data.HFSS_best;
HFSS_Test = HFSS_Test.set_parameter('lline', variables(1), 'mm');
HFSS_Test = HFSS_Test.set_parameter('wpatch2', variables(2), 'mm');
HFSS_Test = HFSS_Test.set_parameter('lpatch2', variables(3), 'mm');
HFSS_Test = HFSS_Test.set_parameter('wpatch3', variables(4), 'mm');
HFSS_Test = HFSS_Test.set_parameter('lpatch3', variables(5), 'mm');
HFSS_Test = HFSS_Test.set_HFSS_parameters('lline');
HFSS_Test = HFSS_Test.set_HFSS_parameters('wpatch2');
HFSS_Test = HFSS_Test.set_HFSS_parameters('lpatch2');
HFSS_Test = HFSS_Test.set_HFSS_parameters('wpatch3');
HFSS_Test = HFSS_Test.set_HFSS_parameters('lpatch3');
HFSS_Test = HFSS_Test.simulate_freq;
HFSS_Test = HFSS_Test.HFSS_export_S_parameters;
HFSS_Test = HFSS_Test.import_S_parameters;
[frequency, S] = HFSS_Test.get_S_param(1, 1);

% freq_pnt_L = find(frequency >= 59e9, 1);
% freq_pnt_U = find(frequency >= 61e9, 1);
cost = max(dBv(S));
if isempty(Opt_Data.cost_min) || cost < Opt_Data.cost_min(end)
    Opt_Data.cost_min(end+1) = cost;
    Opt_Data.HFSS_Best = HFSS_Test;
    Opt_Data.var_best = variables;
else
    Opt_Data.cost_min(end+1) = Opt_Data.cost_min(end);
end
Opt_Data.cost(end+1) = cost;
disp('working');
end

