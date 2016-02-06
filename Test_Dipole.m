clear all
close all


HFSS_Test = HFSS_Tools('Dipole_Opt', 'HFSSDesign1', 1);
HFSS_Test = HFSS_Test.set_solution_setup('Setup1');
HFSS_Test = HFSS_Test.set_solution_sweep_freq('Sweep');
HFSS_Test = HFSS_Test.get_HFSS_parameters;

%%

% HFSS_Test = HFSS_Test.set_parameter('Dp_L', variables(1), 'mm');
% HFSS_Test = HFSS_Test.HFSS_set_parameters('Dp_L');
HFSS_Test = HFSS_Test.simulate_freq;
HFSS_Test = HFSS_Test.get_s_parameters;
HFSS_Test.print_S_param(1);
[frequency, S] = HFSS_Test.get_S_param(1, 1);

freq_pnt_L = find(frequency >= 59e9, 1);
freq_pnt_U = find(frequency >= 61e9, 1);
if isempty(freq_pnt_U)
    freq_pnt_U = length(frequency);
end
cost = mean(dBv(S(freq_pnt_L:freq_pnt_U)));
disp(['Cost : ', num2str(cost)]);