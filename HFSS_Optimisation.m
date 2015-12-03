% clear all;
% close all;

global Opt_Data;
% HFSS_Test = HFSS_Tools('testcell_forLen', 'ltcc_antenna', 1);
% HFSS_Test = HFSS_Test.set_solution_setup('HFSS_Setup_1');
% HFSS_Test = HFSS_Test.set_solution_sweep_freq('Sweep_2');
% HFSS_Test = HFSS_Test.get_HFSS_parameters;
Opt_Data.cost = [];
Opt_Data.cost_min = [];
Opt_Data.mesh_size = [];
Opt_Data.HFSS_best = HFSS_Test;
Opt_Data.var_best = [];

var_start = [0.05,  0.5,    0.9,    1.33,   0.85];
var_min =   [0.05,  0.2,    0.6,    0.8,    0.5];
var_max =   [0.4,   0.8,    1.2,    1.5,    1.2];

PSoptions = psoptimset('PlotFcns', {@opt_plot});
[x, fval] = patternsearch(@run_and_cost_function, var_start, [], [], [], [], var_min, var_max, PSoptions);

Opt_Data.HFSS_best = Opt_Data.HFSS_best.set_HFSS_parameters('All');
