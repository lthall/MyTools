function stop = opt_plot(optimvalues,flag)
%PSPLOTBESTF PlotFcn to plot best function value.
%   STOP = PSPLOTBESTF(OPTIMVALUES,FLAG) where OPTIMVALUES is a structure
%   with the following fields:
%              x: current point X
%      iteration: iteration count
%           fval: function value
%       meshsize: current mesh size
%      funccount: number of function evaluations
%         method: method used in last iteration
%         TolFun: tolerance on function value in last iteration
%           TolX: tolerance on X value in last iteration
%
%   FLAG: Current state in which PlotFcn is called. Possible values are:
%           init: initialization state
%           iter: iteration state
%           done: final state
%
%   STOP: A boolean to stop the algorithm.
%
%   See also PATTERNSEARCH, GA, PSOPTIMSET.


%   Copyright 2003-2005 The MathWorks, Inc.

global Opt_Data;
stop = false;
Opt_Data.mesh_size(end+1) = optimvalues.meshsize;
disp(num2str(Opt_Data.var_best));
switch flag
    case 'init'
        subplot(2,2,1);
            hold on;
            plot(Opt_Data.cost,'b');
            plot(Opt_Data.cost_min,'r');
            xlabel('Simulation','interp','none'); 
            ylabel('Cost','interp','none');
            title('Cost','interp','none');
        subplot(2,2,2);
            plot(Opt_Data.mesh_size,'b');
            xlabel('Simulation','interp','none'); 
            ylabel('Size','interp','none');
            title('Mesh Size','interp','none');
        subplot(2,2,3);
        subplot(2,2,4);
            xlabel('frequency GHz','interp','none'); 
            ylabel('S11 (dB)','interp','none');
            title('S-Parameters','interp','none');
    case 'iter'
        subplot(2,2,1);
            hold on;
            plot(Opt_Data.cost,'b');
            plot(Opt_Data.cost_min,'r');
            xlabel('Simulation','interp','none'); 
            ylabel('Cost','interp','none');
            title('Cost','interp','none');
        subplot(2,2,2);
            plot(Opt_Data.mesh_size,'b');
            xlabel('Simulation','interp','none'); 
            ylabel('Size','interp','none');
            title('Mesh Size','interp','none');
        subplot(2,2,3);
        subplot(2,2,4);
            [frequency, S] = Opt_Data.HFSS_Best.get_S_param(1, 1);
            plot(frequency, 20*log10(abs(S)));
            xlabel('frequency GHz','interp','none'); 
            ylabel('S11 (dB)','interp','none');
            title('S-Parameters','interp','none');
    case 'done'
        subplot(2,2,1);
            hold on;
            plot(Opt_Data.cost,'b');
            plot(Opt_Data.cost_min,'r');
            xlabel('Simulation','interp','none'); 
            ylabel('Cost','interp','none');
            title('Cost','interp','none');
        subplot(2,2,2);
            plot(Opt_Data.mesh_size,'b');
            xlabel('Simulation','interp','none'); 
            ylabel('Size','interp','none');
            title('Mesh Size','interp','none');
        subplot(2,2,3);
        subplot(2,2,4);
            [frequency, S] = Opt_Data.HFSS_Best.get_S_param(1, 1);
            plot(frequency, 20*log10(abs(S)));
            xlabel('frequency GHz','interp','none'); 
            ylabel('S11 (dB)','interp','none');
            title('S-Parameters','interp','none');
end
