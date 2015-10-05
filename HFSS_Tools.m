classdef HFSS_Tools
    %HFSS_TOOLS This provides a range of functions for interfacing with
    %HFSS and analysing the results
    %   Detailed explanation goes here
    
    properties
        project_name;           % HFSS project name
        active_design;          % HFSS design
        source_number;          % Number of sources in HFSS simulation
        source_type;            % Modal, Terminal, Designer
        sources;
        solution_setup;         % solution Setup for fields and S-parameters
        solution_sweep_freq;    % sweep used for S-parameters
        solution_sweep_fields;  % sweep used for fields
        fields_axis = '';       % sweep used for fields
        parameter_list;
        theta_step = 1;         % degrees
        phi_step = 1;           % degrees
        S_matrix;
        beam_patterns;
    end
    
    methods
        %% Setup Methods
        
        function obj = HFSS_Tools(project_name, active_design, source_number)
            %HFSS_TOOLS Initilisation of HFSS_Tools object
            %   Copy text from project_name and active_design HFSS name box
            %   and enter the number of ports in the design.
            
            obj.project_name = project_name;
            obj.active_design = active_design;
            obj.source_number = source_number;
            obj.source_type = 'Modal';
            for aa = 1:obj.source_number
                obj.sources(aa).name = int2str(aa);
                obj.sources(aa).modes = 1;
                obj.sources(aa).excitation_E = 1;
                obj.sources(aa).position = [0,0,0];
                obj.sources(aa).groups = 0;
            end
        end
        
        
        function obj = set_project_name(obj, variable)
            %SET_PROJECT_NAME Set the project name variable
            %    Set the project name variable from the HFSS name box.
            
            obj.project_name = variable;
        end
        
        
        function obj = set_active_design(obj, variable)
            %SET_ACTIVE_DESIGN  Set the design name variable
            %    Set the design name variable from the HFSS name box.
            
            obj.active_design = variable;
        end
        
        
        function obj = set_source_number(obj, variable)
            %SET_SOURCE_NO Set the number of ports
            %    Set the number of sorces variable equal to the number of
            %    ports in HFSS design.
            
            obj.source_number = variable;
        end
        
        
        function obj = set_solution_setup(obj, variable)
            %SET_SOLUTION_SETUP Set the solution setup to use
            %   Both the frequency sweep and fields sweep must use this
            %   solution.
            
            obj.solution_setup = variable;
        end
        
        
        function obj = set_solution_sweep_freq(obj, variable)
            %SET_SOLUTION_SWEEP Set the frequency sweep to use
            %   This is the sweep used to generate s-parameters.
            
            obj.solution_sweep_freq = variable;
        end
        
        
        function obj = set_solution_sweep_fields(obj, variable)
            %SET_SOLUTION_SWEEP_FIELDS Set the fields frequency sweep to use
            %   This is the sweep used to generate the beam patterns.
            
            obj.solution_sweep_fields = variable;
        end
        
        
        function obj = set_fields_axis(obj, variable)
            %SET_SOLUTION_SWEEP_FIELDS Set the fields frequency sweep to use
            %   This is the sweep used to generate the beam patterns.
            
            obj.fields_axis = variable;
        end
        
        
        function obj = set_source_type(obj, variable)
            %SET_SOLUTION_SWEEP_FIELDS Set the fields frequency sweep to use
            %   This is the sweep used to generate the beam patterns.
            
            obj.source_type = variable;
        end
        
        
        function obj = set_angle_resolution(obj, theta_step, phi_step)
            %SET_ANGLE_RESOLUTION Sets the angle resolution for fields
            %   This is the resolution of all beam pattern results.
            %   generated
            
            obj.theta_step = theta_step;
            obj.phi_step = phi_step;
        end
        
        
        function [obj] = set_source_name(obj, port_no, name)
            %SET_SOURCE_EXCITATION Set the excitation of any port
            %   This sets the excitation of port port_no to
            %   excitation_E. port_no should be a positive integer
            %   less than the number of ports. excitation_E is a complex
            %   number representing the magnitude and phase of the
            %   excitation voltage source.
            
            if port_no > length(obj.sources)
                error('That source does not exist');
            else
                obj.sources(port_no).name = name;
            end
        end
        
        
        function [obj] = set_source_groups(obj, port_no, groups)
            %SET_SOURCE_EXCITATION Set the excitation of any port
            %   This sets the excitation of port port_no to
            %   excitation_E. port_no should be a positive integer
            %   less than the number of ports. excitation_E is a complex
            %   number representing the magnitude and phase of the
            %   excitation voltage source.
            
            for aa = port_no
                if port_no > length(obj.sources)
                    error('Source ', int2str(aa), ' does not exist');
                else
                    obj.sources(aa).groups = groups;
                end
            end
        end
        
        
        function [obj] = set_source_excitation_E(obj, port_no, excitation_E)
            %SET_SOURCE_EXCITATION Set the excitation of any port
            %   This sets the excitation of port port_no to
            %   excitation_E. port_no should be a positive integer
            %   less than the number of ports. excitation_E is a complex
            %   number representing the magnitude and phase of the
            %   excitation voltage source.
            
            if port_no > length(obj.sources)
                error('That source does not exist');
            else
                obj.sources(port_no).excitation_E = excitation_E;
            end
        end
        
        
        function [obj] = set_source_position(obj, pnt_antenna, position)
            %SET_SOURCE_POSITION Set the position of each antenna element
            %   This sets the position of port pnt_antenna to
            %   position.  pnt_antenna should be a positive integer
            %   less than the number of ports. position is a three axis
            %   vector representing the position of the antenna in
            %   cartesian coordinates.
            
            if pnt_antenna > length(obj.sources)
                error('That source does not exist');
            else
                obj.sources(pnt_antenna).position = position;
            end
        end
        
        
        function [excitation_E] = get_source_excitation_E(obj, port_no)
            %GET_SOURCE_EXCITATION Get the excitation of any port
            %   This gets the excitation of port port_no. port_no should
            %   be a positive integer less than the number of ports.
            %   excitation_E is a complex number representing the magnitude
            %   and phase of the excitation voltage source.
            
            if port_no > length(obj.sources)
                error('That source does not exist');
            else
                excitation_E = obj.sources(port_no).excitation_E;
            end
        end
        
        
        function [position] = get_source_position(obj, aa)
            %GET_SOURCE_POSITION Get the position of each antenna element
            %   This gets the position of the antenna pnt_antenna.  pnt_antenna
            %   should be a positive integer less than the number of ports.
            %   position is a three axis vector representing the position of
            %   the antenna in cartesian coordinates.
            
            if aa > length(obj.sources)
                error('That source does not exist');
            else
                position = obj.sources(aa).position;
            end
        end
        
        
        %% HFSS Parameter Methods
        
        function obj = get_HFSS_parameters(obj)
            %GET_HFSS_PARAMETERS gets the parameters from the HFSS design
            %   The prompts you to right click on the HFSS parameter window
            %   and select copy to clipboard. Matlab then populates the
            %   paramers from this clipboard data.
            
            if isempty(obj.project_name) || isempty(obj.active_design)
                error('You must supply all HFSS information');
            else
                hfss_script = [];
                [ hfss_script ] = obj.HFSS_Open( hfss_script );
                fileID = fopen('run_me_get_param.vbs', 'w');
                fprintf(fileID, hfss_script);
                fclose(fileID);
                system('run_me_get_param.vbs');
            end
            
            questdlg('Right click HFSS parameters and copy to clipboard. Then press OK.', 'Click when ready', 'OK', 'OK');
            parameter_string = clipboard('paste');
            param_data = textscan(parameter_string, '%s', 'delimiter','\n');
            tab = sprintf('\t');
            obj.parameter_list = [];
            for aa = 2:length(param_data{1})
                line_data = param_data{1}(aa);
                tab_pnt = [];
                for bb = 1:length(line_data{1})
                    if strcmp(tab,line_data{1}(bb))
                        tab_pnt(end+1) = bb;
                    end
                end
                obj.parameter_list(aa-1).name = line_data{1}(1:tab_pnt(1)-1);
                
                if strcmp(line_data{1}(tab_pnt(1)+1), '"')
                    obj.parameter_list(aa-1).type = 'eqn';
                else
                    obj.parameter_list(aa-1).type = 'var';
                end
                
                line_pnt_start_unit = tab_pnt(3)+1;
                while ~isletter(line_data{1}(line_pnt_start_unit))
                    line_pnt_start_unit = line_pnt_start_unit + 1;
                end
                value_read = textscan(line_data{1}(tab_pnt(3)+1:line_pnt_start_unit-1), '%f');
                obj.parameter_list(aa-1).value = abs(value_read{1});
                if tab_pnt(4) > line_pnt_start_unit
                    value_read = textscan(line_data{1}(line_pnt_start_unit:tab_pnt(4)-1), '%s');
                    obj.parameter_list(aa-1).units = value_read{1}{1};
                else
                    obj.parameter_list(aa-1).units = 'none';
                end
            end
        end
        
        
        function [ obj ] = set_HFSS_parameters(obj, name)
            %SET_HFSS_PARAMETERS Sets the value of a parameter in HFSS
            %   This writes the value of the parameter specified by the
            %   input variable name. By setting name to 'All' matlab
            %   updates all parameter values in HFSS to be the same as
            %   stored in Matlab.
            
            if isempty(obj.project_name) || isempty(obj.active_design)
                error('You must supply all HFSS information');
            else
                hfss_script = [];
                [ hfss_script ] = obj.HFSS_Open( hfss_script );
                if strcmp('All', name)
                    for pnt = 1:length(obj.parameter_list)
                        if strcmp('var',obj.parameter_list(pnt).type)
                            [ hfss_script ] = HFSS_Save_Parameters(obj, hfss_script, obj.parameter_list(pnt).name, obj.parameter_list(pnt).value, obj.parameter_list(pnt).units );
                        end
                    end
                else
                    pnt = 0;
                    for aa = 1:length(obj.parameter_list)
                        if strcmp(obj.parameter_list(aa).name, name)
                            pnt = aa;
                        end
                    end
                    if pnt == 0
                        error('parameter does not exist');
                    else
                        [ hfss_script ] = obj.HFSS_Save_Parameters(hfss_script, obj.parameter_list(pnt).name, obj.parameter_list(pnt).value, obj.parameter_list(pnt).units );
                    end
                end
                fileID = fopen('run_me_set_param.vbs', 'w');
                fprintf(fileID, hfss_script);
                fclose(fileID);
                system('run_me_set_param.vbs');
            end
        end
        
        
        function [value] = get_parameter(obj, name)
            %GET_PARAMETER returns the parameter value
            %   This returns the parameter value of the parameter specified
            %   by the name input. All values are returned in base SI
            %   units, radians and meters.
            
            pnt = 0;
            for aa = 1:length(obj.parameter_list)
                if strcmp(obj.parameter_list(aa).name, name)
                    pnt = aa;
                end
            end
            if pnt > 0
                switch obj.parameter_list(pnt).units
                    case 'm'
                        unit = 1;
                    case {'mm'}
                        unit = 1e-3;
                    case {'um'}
                        unit = 1e-6;
                    case {'u'}
                        unit = 1e-6;
                    case {'in'}
                        unit = inch;
                    otherwise
                        unit = 1;
                end
                value = obj.parameter_list(pnt).value*unit;
            else
                error('Parameter does not exist')
            end
        end
        
        
        function [obj] = set_parameter(obj, name, value, units)
            %SET_PARAMETER sets the value of a parameter
            %   This function sets the value of a parameter in matlab using
            %   standard units. To write these values to HFSS you must use
            %   the [ obj ] = set_HFSS_parameters(obj, name) function.
            
            pnt = 0;
            for aa = 1:length(obj.parameter_list)
                if strcmp(obj.parameter_list(aa).name, name)
                    pnt = aa;
                end
            end
            if pnt > 0
                obj.parameter_list(pnt).value = value;
                obj.parameter_list(pnt).units = units;
            else
                error('Parameter does not exist')
            end
        end
        
        
        function list_parameters(obj)
            %LIST_PARAMETERS outputs a list of all variables
            
            for aa = 1:length(obj.parameter_list)
                fprintf([obj.parameter_list(aa).name, '\t Value: ', num2str(obj.parameter_list(aa).value), '\t Units: ', obj.parameter_list(aa).units, '\t Type: ', obj.parameter_list(aa).type, '\n']);
            end
        end
        
        
        %% HFSS Simulation Methods
        
        function obj = simulate_freq(obj)
            %SIMULATE_FREQ starts a frequency simulation
            %   This function starts the sweep specified to be used to
            %   generate s-parameters. If you have not specified a
            %   frequency sweep it will use the single point defined by
            %   your solution setup.
            
            if isempty(obj.project_name) || isempty(obj.active_design) || isempty(obj.solution_setup) || isempty(obj.source_number)
                error('You must supply all HFSS information');
            else
                hfss_script = [];
                [ hfss_script ] = obj.HFSS_Open( hfss_script );
                [ hfss_script ] = obj.HFSS_Simulate_freq( hfss_script );
                fileID = fopen('run_me_simulate.vbs', 'w');
                fprintf(fileID, hfss_script);
                fclose(fileID);
                system('run_me_simulate.vbs');
            end
        end
        
        
        function obj = simulate_fields(obj)
            %SIMULATE_FIELDS starts a fields frequency simulation
            %   This function starts the sweep specified to be used to
            %   generate fields reports. If you have not specified a
            %   frequency sweep it will use the single point defined by
            %   your solution setup.
            
            if isempty(obj.project_name) || isempty(obj.active_design) || isempty(obj.solution_setup) || isempty(obj.source_number)
                error('You must supply all HFSS information');
            else
                hfss_script = [];
                [ hfss_script ] = obj.HFSS_Open( hfss_script );
                [ hfss_script ] = obj.HFSS_Simulate_fields( hfss_script );
                fileID = fopen('run_me_simulate.vbs', 'w');
                fprintf(fileID, hfss_script);
                fclose(fileID);
                system('run_me_simulate.vbs');
            end
        end
        
        
        %% HFSS Export and Import Result Methods
        
        function obj = HFSS_export_S_parameters(obj)
            %HFSS_EXPORT_S_PARAMETERS Export s-parameters in HFSS to file
            %   This function exports the s-parameters associated with your
            %   frequency sweep and exports them to a sNp file in your
            %   working directory.
            
            if isempty(obj.project_name) || isempty(obj.active_design) || isempty(obj.solution_setup) || isempty(obj.source_number)
                error('You must supply all HFSS information');
            else
                hfss_script = [];
                [ hfss_script ] = obj.HFSS_Open( hfss_script );
                file_name = [pwd,'\sparam_', obj.project_name, '_', obj.active_design, '.', 's', int2str(obj.source_number), 'p'];
                [ hfss_script ] = obj.HFSS_S_Parameters( hfss_script, file_name );
                fileID = fopen('run_me_s_param.vbs', 'w');
                fprintf(fileID, hfss_script);
                fclose(fileID);
                system('run_me_s_param.vbs');
            end
        end
        
        
        function obj = DESIGNER_export_S_parameters(obj)
            %HFSS_EXPORT_S_PARAMETERS Export s-parameters in HFSS to file
            %   This function exports the s-parameters associated with your
            %   frequency sweep and exports them to a sNp file in your
            %   working directory.
            
            questdlg(['Please export S-Parameters to : ', pwd,'\sparam.', 's', int2str(obj.source_number), 'p'], 'Click when ready', 'OK', 'OK');
        end
        
        
        function obj = import_S_parameters(obj, file_name)
            %IMPORT_S_PARAMETERS Import s-parameters from SnP file
            %   This function imports the SnP file specified by the file
            %   name file_name. If you don't specify a file name the code
            %   assumes the SnP file has been exported from HFSS using this
            %   toolbox.
            
            if isempty(obj.project_name) || isempty(obj.active_design) || isempty(obj.solution_setup) || isempty(obj.source_number)
                error('You must supply all HFSS information');
            else
                if nargin == 2
                    file_name_SnP = file_name;
                else
                    file_name_SnP = [pwd,'\sparam_', obj.project_name, '_', obj.active_design, '.', 's', int2str(obj.source_number), 'p'];
                end
                obj.S_matrix = obj.import_SnP(file_name_SnP);
            end
        end
        
        
        function obj = HFSS_export_fields_data(obj)
            %HFSS_EXPORT_FIELDS_DATA Export fields data in HFSS to file
            %   This function exports the fields data associated with your
            %   fields frequency sweep. There is one set of field data
            %   exported per port in HFSS. Each set of field data is
            %   exported to your working directory as a .csv file.
            
            if isempty(obj.project_name) || isempty(obj.active_design) || isempty(obj.solution_setup) || isempty(obj.source_number)
                error('You must supply all HFSS information');
            else
                hfss_script = [];
                [ hfss_script ] = obj.HFSS_Open( hfss_script );
                switch obj.source_type
                    case 'Modal'
                        [ hfss_script ] = obj.HFSS_Set_Mode_Sources( hfss_script );
                        [ hfss_script ] = obj.HFSS_Add_Rad_Sphere( hfss_script );
                        [ hfss_script ] = obj.HFSS_Add_Rad_Report( hfss_script );
                    case 'Terminal'
                        [ hfss_script ] = obj.HFSS_Set_Terminal_Sources( hfss_script );
                        [ hfss_script ] = obj.HFSS_Add_Rad_Sphere( hfss_script );
                        [ hfss_script ] = obj.HFSS_Add_Rad_Report( hfss_script );
                    case 'Designer'
                        [ hfss_script ] = obj.DESIGNER_Set_Sources( hfss_script );
                        [ hfss_script ] = obj.DESIGNER_Add_Rad_Report( hfss_script );
                end
                for aa = 1:obj.source_number
                    obj.sources(aa).excitation_E = 1;
                    for bb = 1:obj.source_number
                        if bb ~= aa
                            obj.sources(bb).excitation_E = 0;
                        end
                    end
                    file_name = [pwd,'\fields_', obj.project_name, '_', obj.active_design, '_', int2str(aa)];
                    switch obj.source_type
                        case 'Modal'
                            [ hfss_script ] = obj.HFSS_Set_Mode_Sources( hfss_script );
                        case 'Terminal'
                            [ hfss_script ] = obj.HFSS_Set_Terminal_Sources( hfss_script );
                        case 'Designer'
                            [ hfss_script ] = obj.DESIGNER_Set_Sources( hfss_script );
                    end
                    [ hfss_script ] = obj.HFSS_Save_Rad_Report( hfss_script, file_name );
                end
                switch obj.source_type
                    case 'Modal'
                        [ hfss_script ] = obj.HFSS_Remove_Matlab( hfss_script );
                    case 'Terminal'
                        [ hfss_script ] = obj.HFSS_Remove_Matlab( hfss_script );
                    case 'Designer'
                        [ hfss_script ] = obj.DESIGNER_Remove_Matlab( hfss_script );
                        
                end
                fileID = fopen('run_me_get_fields.vbs', 'w');
                fprintf(fileID, hfss_script);
                fclose(fileID);
                system('run_me_get_fields.vbs');
            end
        end
        
        
        function [ obj ] = import_fields_data(obj, file_name)
            %IMPORT_FIELDS_DATA Import fields data from csv files
            %   This function imports the csv files specified by the file
            %   name file_name. If you don't specify a file name the code
            %   assumes the csv files have been exported from HFSS using this
            %   toolbox.
            
            for aa = 1:obj.source_number
                if nargin == 2
                    file_name_fields = file_name;
                else
                    file_name_fields = [pwd,'\fields_', obj.project_name, '_', obj.active_design, '_', int2str(aa)];
                end
                antennas(aa).parameters = obj.import_fields(file_name_fields);
            end
            
            pnt_freq = 0;
            for aa = 1:length(antennas(1).parameters(1).variable)
                if strcmp(antennas(1).parameters(1).variable(aa).name, 'Freq') || strcmp(antennas(1).parameters(1).variable(aa).name, 'F')
                    pnt_freq = aa;
                    break;
                end
            end
            
            beams(1).frequency = antennas(1).parameters(1).variable(pnt_freq).value;
            freq(1).pnts = [];
            found = 0;
            for aa = 1:length(antennas(1).parameters)
                for bb = 1:length(beams)
                    if beams(bb).frequency == antennas(1).parameters(aa).variable(pnt_freq).value
                        found = 1;
                        freq(bb).pnts(end+1) = aa;
                    end
                end
                if found == 0
                    beams(end+1).frequency = antennas(1).parameters(aa).variable(pnt_freq).value;
                    freq(end+1).pnts(1) = aa;
                end
                found = 0;
            end
            
            for aa = 1:obj.source_number
                for bb = 1:length(antennas(aa).parameters)
                    antennas(aa).phi(bb,:) = ones(size(antennas(aa).parameters(bb).output(1).value)) * antennas(aa).parameters(bb).variable(2).value;
                    antennas(aa).theta(bb,:) = antennas(aa).parameters(bb).output(1).value;
                    antennas(aa).rE_phi(bb,:) = antennas(aa).parameters(bb).output(2).value;
                    antennas(aa).rE_theta(bb,:) = antennas(aa).parameters(bb).output(3).value;
                end
            end
            
            phi = antennas(1).phi(freq(1).pnts,:);
            theta = antennas(1).theta(freq(1).pnts,:);
            for aa = 1:length(beams)
                beams(aa).phi = phi;
                beams(aa).theta = theta;
                for bb = 1:length(antennas)
                    beams(aa).antenna(bb).rE_phi = antennas(bb).rE_phi(freq(aa).pnts,:);
                    beams(aa).antenna(bb).rE_theta = antennas(bb).rE_theta(freq(aa).pnts,:);
                end
            end
            obj.beam_patterns = beams;
        end
        
        
        %% Access Methods
        
        function [frequency, S, legend_text] = print_S_param(obj, fig_no)
            %PRINT_S_PARAM displays the S-parameter data
            %   This function displays and returns the S-parameter data.
            
            S_pnt = 0;
            for a = 1:obj.source_number
                for b = 1:obj.source_number
                    S_pnt = S_pnt + 1;
                    for c = 1:length(obj.S_matrix)
                        frequency(c) = obj.S_matrix(c).frequency;
                        S(c,S_pnt) = obj.S_matrix(c).S(a,b);
                        if a == b
                            Snn(c,a) = obj.S_matrix(c).S(a,b);
                            Snn_legend_text{a} = ['S_{', int2str(a), int2str(b), '}'];
                        end
                    end
                    legend_text{S_pnt} = ['S_{', int2str(a), int2str(b), '}'];
                end
            end
            
            if ~isempty(fig_no)
                figure(fig_no)
                subplot(2,1,1);
                hold on
                plot(frequency/10^9,20*log10(abs(S)));
                ylabel('Mag dB');
                grid on;
                subplot(2,1,2);
                hold on
                plot(frequency/10^9,angle(S)/HFSS_Tools.degrees);
                ylabel('Angle Deg');
                grid on;
                xlabel('frequency GHz');
                legend(legend_text);
                
                figure(fig_no+1)
                clf
                obj.MrSmith;
                plot(Snn);
                
                figure(fig_no+2)
                clf
                subplot(2,1,1);
                hold on
                plot(frequency/10^9,20*log10(abs(Snn)));
                ylabel('Mag dB');
                grid on;
                subplot(2,1,2);
                hold on
                plot(frequency/10^9,angle(Snn)/HFSS_Tools.degrees);
                ylabel('Angle Deg');
                grid on;
                xlabel('frequency GHz');
                legend(Snn_legend_text);
            end
        end
        
        
        function [frequency, S, legend_text] = print_S_param_all(obj)
            %PRINT_S_PARAM displays the S-parameter data
            %   This function displays and returns the S-parameter data.
            
            S_pnt = 0;
            for a = 1:obj.source_number
                for b = 1:obj.source_number
                    S_pnt = S_pnt + 1;
                    for c = 1:length(obj.S_matrix)
                        frequency(c) = obj.S_matrix(c).frequency;
                        S(c,S_pnt) = obj.S_matrix(c).S(a,b);
                        if a == b
                            Snn(c,a) = obj.S_matrix(c).S(a,b);
                            Snn_legend_text{a} = ['S_{', int2str(a), int2str(b), '}'];
                        end
                    end
                    legend_text{S_pnt} = ['S_{', int2str(a), int2str(b), '}'];
                end
            end
            
            figure
            subplot(2,1,1);
            hold on
            plot(frequency/10^9,20*log10(abs(S)));
            ylabel('Mag dB');
            grid on;
            subplot(2,1,2);
            hold on
            plot(frequency/10^9,angle(S)/HFSS_Tools.degrees);
            ylabel('Angle Deg');
            grid on;
            xlabel('frequency GHz');
            legend(legend_text);
            
            figure
            obj.MrSmith;
            plot(Snn);
            
            figure
            subplot(2,1,1);
            hold on
            plot(frequency/10^9,20*log10(abs(Snn)));
            ylabel('Mag dB');
            grid on;
            subplot(2,1,2);
            hold on
            plot(frequency/10^9,angle(Snn)/HFSS_Tools.degrees);
            ylabel('Angle Deg');
            grid on;
            xlabel('frequency GHz');
            legend(Snn_legend_text);
        end
        
        
        function [frequency, S, legend_text] = print_S_param_group(obj, group)
            %PRINT_S_PARAM displays the S-parameter data
            %   This function displays and returns the S-parameter data.
            
            S_pnt = 0;
            legend_pnt = 0;
            for a = 1:obj.source_number
                for b = 1:obj.source_number
                    if (ismember(group, obj.sources(a).groups) && ismember(group, obj.sources(b).groups))
                        S_pnt = S_pnt + 1;
                        for c = 1:length(obj.S_matrix)
                            frequency(c) = obj.S_matrix(c).frequency;
                            S(c,S_pnt) = obj.S_matrix(c).S(a,b);
                            if a == b
                                if c == 1
                                    legend_pnt = legend_pnt + 1;
                                end
                                Snn(c,a) = obj.S_matrix(c).S(a,b);
                                Snn_legend_text{legend_pnt} = ['S_{', int2str(a), int2str(b), '}'];
                            end
                        end
                        legend_text{S_pnt} = ['S_{', int2str(a), int2str(b), '}'];
                    end
                end
            end
            
            figure
            subplot(2,1,1);
            hold on
            plot(frequency/10^9,20*log10(abs(S)));
            ylabel('Mag dB');
            grid on;
            subplot(2,1,2);
            hold on
            plot(frequency/10^9,angle(S)/HFSS_Tools.degrees);
            ylabel('Angle Deg');
            grid on;
            xlabel('frequency GHz');
            legend(legend_text);
            
            figure
            obj.MrSmith;
            plot(Snn);
            
            figure
            subplot(2,1,1);
            hold on
            plot(frequency/10^9,20*log10(abs(Snn)));
            ylabel('Mag dB');
            grid on;
            subplot(2,1,2);
            hold on
            plot(frequency/10^9,angle(Snn)/HFSS_Tools.degrees);
            ylabel('Angle Deg');
            grid on;
            xlabel('frequency GHz');
            legend(Snn_legend_text);
        end
        
        
        function [frequency, S] = get_S_param(obj, nn, mm)
            %GET_S_PARAM Returns S-parameter data
            %   Returns the S-parameter data specified by the variables nn
            %   and mm and all frequency points.
            
            for aa = 1:length(obj.S_matrix)
                frequency(aa) = obj.S_matrix(aa).frequency;
                S(aa) = obj.S_matrix(aa).S(nn,mm);
            end
        end
        
        
        function [frequency] = get_S_freq(obj)
            %GET_S_FREQ Returns the vector of frequency points used by the
            %S-parameter data
            
            frequency = ones(length(obj.S_matrix),1);
            for aa = 1:length(obj.S_matrix)
                frequency(aa) = obj.S_matrix(aa).frequency;
            end
        end
        
        
        function [frequency] = get_fields_freq(obj)
            %GET_FIELDS_FREQ Returns the vector of frequency points used by the
            %fields data
            
            frequency = ones(length(obj.beam_patterns),1);
            for aa = 1:length(obj.beam_patterns)
                frequency(aa) = obj.beam_patterns(aa).frequency;
            end
        end
        
        
        %% Processing Methods
        
        function [ return_voltage_wave, frequency ] = array_return_voltage_wave( obj, group )
            %ARRAY_RETURN_VOLTAGE_WAVE returns the return voltage wave
            %for an array excitation
            %   The return voltage wave is the return wave present at each
            %   port of the array. The value of these parameters depends on the
            %   excitation of each element of the array.
            
            frequency = ones(length(obj.S_matrix),1);
            return_voltage_wave = ones(length(obj.S_matrix),obj.source_number);
            for aa = 1:length(obj.S_matrix);
                source_E = ones(obj.source_number,1);
                for bb = 1:obj.source_number;
                    if isempty(group) || ismember(group, obj.sources(bb).groups)
                        if length(obj.S_matrix) == length(obj.sources(bb).excitation_E)
                            source_E(bb) = obj.sources(bb).excitation_E(aa);
                        elseif length(obj.sources(bb).excitation_E) == 1
                            source_E(bb) = obj.sources(bb).excitation_E(1);
                        else
                            disp('You have not provided an exciation at each frequency point, using closest point interpolation');
                            if obj.S_matrix(aa).frequency <= obj.sources(bb).frequency(1)
                                source_E(bb) = obj.sources(bb).excitation_E(1);
                            elseif obj.S_matrix(aa).frequency >= obj.sources(bb).frequency(end)
                                source_E(bb) = obj.sources(bb).excitation_E(end);
                            else
                                source_E(bb) = interp1(obj.sources(bb).frequency, obj.sources(bb).excitation_E, obj.S_matrix(aa).frequency, 'nearest');
                            end
                        end
                    else
                        source_E(bb) = 0;
                    end
                end
                frequency(aa) = obj.S_matrix(aa).frequency;
                return_voltage_wave(aa,:) = obj.S_matrix(aa).S * source_E;
            end
        end
        
        
        function [ return_power, frequency ] = array_return_power( obj, group )
            %ARRAY_RETURN_POWER returns the total return power of the array
            %   The array return power is the total power reflected back
            %   into the ports of the array. The return power depends on the
            %   excitation of each element of the array.
            
            frequency = ones(length(obj.S_matrix),1);
            return_voltage_wave = ones(length(obj.S_matrix),obj.source_number);
            return_power = ones(length(obj.S_matrix),1);
            for aa = 1:length(obj.S_matrix);
                source_E = ones(obj.source_number,1);
                for bb = 1:obj.source_number;
                    if isempty(group) || ismember(group, obj.sources(bb).groups)
                        if length(obj.S_matrix) == length(obj.sources(bb).excitation_E)
                            source_E(bb) = obj.sources(bb).excitation_E(aa);
                        elseif length(obj.sources(bb).excitation_E) == 1
                            source_E(bb) = obj.sources(bb).excitation_E(1);
                        else
                            disp('You have not provided an exciation at each frequency point, using closest point interpolation');
                            if obj.S_matrix(aa).frequency <= obj.sources(bb).frequency(1)
                                source_E(bb) = obj.sources(bb).excitation_E(1);
                            elseif obj.S_matrix(aa).frequency >= obj.sources(bb).frequency(end)
                                source_E(bb) = obj.sources(bb).excitation_E(end);
                            else
                                source_E(bb) = interp1(obj.sources(bb).frequency, obj.sources(bb).excitation_E, obj.S_matrix(aa).frequency, 'nearest');
                            end
                        end
                    else
                        source_E(bb) = 0;
                    end
                end
                frequency(aa) = obj.S_matrix(aa).frequency;
                return_voltage_wave(aa,:) = obj.S_matrix(aa).S * source_E';
                return_power(aa) = sum(abs(return_voltage_wave(aa,:)).^2)/sum(abs(source_E).^2);
            end
        end
        
        
        function [reflection_coefficient, frequency ] = array_active_Snn( obj, group )
            %ARRAY_ACTIVE_SNN returns the active s-parameters for a given
            %array excitation
            %   The active s-parameters represent the apparent impedance
            %   seen at each port of the array when the array is being
            %   excited. The value of these parameters depends on the
            %   excitation of each element of the array.
            
            frequency = ones(length(obj.S_matrix),1);
            return_voltage_wave = ones(length(obj.S_matrix),obj.source_number);
            reflection_coefficient = ones(length(obj.S_matrix),obj.source_number);
            for aa = 1:length(obj.S_matrix);
                source_E = ones(obj.source_number,1);
                for bb = 1:obj.source_number;
                    if isempty(group) || ismember(group, obj.sources(bb).groups)
                        if length(obj.S_matrix) == length(obj.sources(bb).excitation_E)
                            source_E(bb) = obj.sources(bb).excitation_E(aa);
                        elseif length(obj.sources(bb).excitation_E) == 1
                            source_E(bb) = obj.sources(bb).excitation_E(1);
                        else
                            disp('You have not provided an exciation at each frequency point, using closest point interpolation');
                            if obj.S_matrix(aa).frequency <= obj.sources(bb).frequency(1)
                                source_E(bb) = obj.sources(bb).excitation_E(1);
                            elseif obj.S_matrix(aa).frequency >= obj.sources(bb).frequency(end)
                                source_E(bb) = obj.sources(bb).excitation_E(end);
                            else
                                source_E(bb) = interp1(obj.sources(bb).frequency, obj.sources(bb).excitation_E, obj.S_matrix(aa).frequency, 'nearest');
                            end
                        end
                    else
                        source_E(bb) = 0;
                    end
                end
                frequency(aa) = obj.S_matrix(aa).frequency;
                return_voltage_wave(aa,:) = obj.S_matrix(aa).S * source_E;
                reflection_coefficient(aa,:) = return_voltage_wave(aa,:)./source_E.';
            end
        end
        
        
        function [ obj ] = set_source_excitation_direction( obj, theta, phi, frequency )
            %SET_SOURCE_EXCITATION_DIRECTION Set sources based on a
            %direction vector
            %   This function sets the phase of each source assuming
            %   identical elements and the position specified for each
            %   port. The phase is set based on the frequency vector. If
            %   only on frequency point is specified it is assuemd the
            %   phase is set at this point and does not vary over
            %   frequency.
            
            direction = obj.direction( theta, phi );
            for aa = 1:obj.source_number;
                for bb = 1:length(frequency)
                    wave_length = obj.c0/frequency(bb);
                    obj.sources(aa).frequency(bb) = frequency(bb);
                    obj.sources(aa).excitation_E(bb) = exp(1j*((-2*pi* dot(obj.sources(aa).position, direction)/wave_length)));
                end
            end
        end
        
        
        function [ obj ] = set_source_max_rE_theta( obj, theta, phi )
            %SET_SOURCE_MAX_RE_THETA Set sources to acheve maximum gain in the
            %direction specified by theta phi
            %   This function sets the phase of each source based on the
            %   beam pattern of each element at each frequency point. The
            %   magnitude of each excitation is set to 1 and the phase the
            %   phase of the complex congergate rE_theta at the desired angle.
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            for aa = 1:length(obj.beam_patterns)
                for bb = 1:obj.source_number
                    X = obj.beam_patterns(aa).theta;
                    Y = obj.beam_patterns(aa).phi;
                    Z = obj.beam_patterns(aa).antenna(bb).rE_theta;
                    ZI = interp2(X,Y,Z,theta,phi);
                    obj.sources(bb).frequency(aa) = obj.beam_patterns(aa).frequency;
                    obj.sources(bb).excitation_E(aa) = (ZI/abs(ZI))';
                end
            end
        end
        
        
        function [ obj ] = set_source_max_rE_phi( obj, theta, phi )
            %SET_SOURCE_MAX_RE_PHI Set sources to acheve maximum gain in the
            %direction specified by theta phi
            %   This function sets the phase of each source based on the
            %   beam pattern of each element at each frequency point. The
            %   magnitude of each excitation is set to 1 and the phase the
            %   phase of the complex congergate rE_theta at the desired angle.
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            for aa = 1:length(obj.beam_patterns)
                for bb = 1:obj.source_number
                    X = obj.beam_patterns(aa).theta;
                    Y = obj.beam_patterns(aa).phi;
                    Z = obj.beam_patterns(aa).antenna(bb).rE_phi;
                    ZI = interp2(X,Y,Z,theta,phi);
                    obj.sources(bb).frequency(aa) = obj.beam_patterns(aa).frequency;
                    obj.sources(bb).excitation_E(aa) = (ZI/abs(ZI))';
                end
            end
        end
        
        
        function [ obj ] = set_source_max_rE( obj, theta, phi, polarization )
            
            if strcmp(polarization,'theta')
                obj = obj.set_source_max_rE_theta( theta, phi );
            elseif strcmp(polarization,'phi')
                obj = obj.set_source_max_rE_phi( theta, phi );
            end
            
        end
        
        
        function [ beam_pattern ] = array_beam_pattern( obj, group )
            %ARRAY_BEAM_PATTERN generates the array beam pattern
            %   This function calculates the array beam pattern based on
            %   excitation_E. If excitation_E is defined at only one
            %   frequency point then the same value is used by all
            %   frequency points.
            
            for aa = 1:length(obj.beam_patterns)
                rE_phi = zeros(size(obj.beam_patterns(aa).phi));
                rE_theta = zeros(size(obj.beam_patterns(aa).theta));
                tx_power = 0;
                for bb = 1:obj.source_number
                    if ismember(group, obj.sources(bb).groups)
                        if length(obj.beam_patterns) == length(obj.sources(bb).excitation_E)
                            source_E = obj.sources(bb).excitation_E(aa);
                        elseif length(obj.sources(bb).excitation_E) == 1
                            source_E = obj.sources(bb).excitation_E(1);
                        else
                            error('You must define an excitation for each frequency point in the s-parameters or just one excitation');
                        end
                        rE_phi = rE_phi + obj.beam_patterns(aa).antenna(bb).rE_phi * source_E;
                        rE_theta = rE_theta + obj.beam_patterns(aa).antenna(bb).rE_theta * source_E;
                        tx_power = tx_power + abs(obj.sources(bb).excitation_E(aa)^2);
                    end
                end
                beam_pattern(aa).frequency = obj.beam_patterns(aa).frequency;
                beam_pattern(aa).phi = obj.beam_patterns(aa).phi;
                beam_pattern(aa).theta = obj.beam_patterns(aa).theta;
                beam_pattern(aa).rE_phi = rE_phi;
                beam_pattern(aa).rE_theta = rE_theta;
                beam_pattern(aa).realized_gain_phi = (2*pi*abs(rE_phi.^2))./(obj.wave_impedance()*tx_power);
                beam_pattern(aa).realized_gain_theta = (2*pi*abs(rE_theta.^2))./(obj.wave_impedance()*tx_power);
                beam_pattern(aa).realized_gain_total = (2*pi*abs(rE_phi.^2+rE_theta.^2))./(obj.wave_impedance()*tx_power);
            end
        end
        
    end
    
    
    methods(Static = true)
        %% File Import Methods
        
        function [ direction ] = direction( theta, phi )
            %DIRECTION calculates a direction vector based on theta and phi
            %   direction is a three element unit vector in cartisian
            %   corrdinates.
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            direction = [cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)];
        end
        
        
        function [S_matrix] = import_SnP(file_name)
            %IMPORT_SnP Imports the S_matrix data
            % This function imports any touchstone SnP file.
            
            Digits = 0;
            while ~isempty(str2num(file_name(length(file_name)-1-Digits)))
                Digits = Digits +1;
            end
            n = str2num(file_name(length(file_name)-(Digits:-1:1)));
            
            FID = fopen(file_name);
            if FID == -1
                error('You got the file name wrong');
                S_matrix = -1;
                return
            end
            
            fseek(FID, 0, 'eof');
            Data_End = ftell(FID);
            frewind(FID);
            
            input_Line = fgetl(FID);
            while ischar(input_Line) && ~strncmp(input_Line, {'#'}, 1)
                input_Line = fgetl(FID);
                if isempty(input_Line)
                    input_Line = ' ';
                end
            end
            if input_Line == -1
                error('Something is wrong with the file')
                S_matrix = -1;
                return
            end
            
            Operators = textscan(input_Line,'%q','delimiter',' ');
            
            switch Operators{1}{2}
                case {'Hz', 'hz', 'HZ'}
                    frequency_Multiplier = 1;
                case {'GHz', 'GHZ'}
                    frequency_Multiplier = 1e9;
                otherwise
                    error('Something is wrong with the file')
                    S_matrix = -1;
                    delete(h);
                    return
            end
            
            switch Operators{1}{3}
                case {'S'}
                    Perameter_Type = 1; % S Perameters
                otherwise
                    error('Something is wrong with the file')
                    S_matrix = -1;
                    delete(h);
                    return
            end
            
            switch Operators{1}{4}
                case {'MA', 'ma'}
                    Magnitude_Scaler = 0; % Magnitude
                case {'DB', 'dB', 'db'}
                    Magnitude_Scaler = 1; % DB
                case {'RI', 'ri', 'Ri', 'rI' }
                    Magnitude_Scaler = 2; % RI
                otherwise
                    error('Something is wrong with the file')
                    S_matrix = -1;
                    delete(h);
                    return
            end
            
            Angle_Multiplier = HFSS_Tools.degrees;
            
            input_Line = fgetl(FID);
            if isempty(input_Line)
                input_Line = '!';
            end
            count_F = 0;
            count_S = 0;
            Time_Loop = [];
            Data_Start = ftell(FID);
            tic
            Next_Update = 1;
            while input_Line ~= -1
                if input_Line(1) ~= '!'
                    [Data] = textscan(input_Line,'%n');
                    Data = Data{1};
                    if ~isempty(Data);
                        if count_S == 0
                            count_F = count_F + 1;
                            S_matrix(count_F).frequency = Data(1)*frequency_Multiplier;
                            Data = Data(2:length(Data));
                        end
                        for a = 1:2:length(Data)
                            count_S = count_S + 1;
                            switch Magnitude_Scaler
                                case 0
                                    S{count_F}(count_S) = Data(a)*exp(1j*Data(a+1)*Angle_Multiplier);
                                case 1
                                    S{count_F}(count_S) = 10^(Data(a)/20)*exp(1j*Data(a+1)*Angle_Multiplier);
                                case 2
                                    S{count_F}(count_S) = Data(a)+1j*Data(a+1);
                            end
                        end
                    end
                end
                if count_S == n^2
                    count_S = 0;
                end
                input_Line = fgetl(FID);
                if isempty(input_Line)
                    input_Line = '!';
                end
                if rem(count_F,10) == 0;
                    Next_Update = min([Next_Update*2,2^7]);
                    Time_Loop(end+1) = toc;
                    tic
                    Time_Left = (Data_End-ftell(FID)) * sum(Time_Loop)/(ftell(FID)-Data_Start);
                end
            end
            
            for a = 1:length(S)
                count_S = 0;
                for b = 1:n
                    for c = 1:n
                        count_S = count_S + 1;
                        S_matrix(a).S(b,c) = S{a}(count_S);
                    end
                end
            end
            fclose(FID);
        end
        
        
        function [parameters] = import_fields(file_name)
            %IMPORT_FIELDS Imports the fields data from HFSS csv file
            % This function imports HFSS table exported as a .csv file.
            
            FID = fopen([file_name, '.csv']);
            if FID == -1
                error('You got the file name wrong');
            end
            
            input_Line = fgetl(FID);
            Outputs_And_Variables = textscan(input_Line,'%q','delimiter',',');
            
            for aa = 1:length(Outputs_And_Variables{1})
                output = textscan(Outputs_And_Variables{1}{aa},'%q','delimiter','-');
                for bb = 3:length(output{1})
                    output{1}{2} = [output{1}{2},'-',output{1}{bb}];
                end
                output_temp = textscan(output{1}{1},'%q','delimiter',' ');
                results(aa).name = output_temp{1};
                if length(output_temp{1}) == 2
                    output_temp{1}{2} = output_temp{1}{2}(2:length(output_temp{1}{2})-1);
                    unit = output_temp{1}{2};
                else
                    unit = '';
                end
                results(aa).units = unit_scailer(unit);
                
                if length(output{1})>1
                    output_temp = textscan(output{1}{2},'%q','delimiter',', ');
                    for bb = 1:length(output_temp{1})
                        text = output_temp{1}{bb};
                        if ~isempty(text)
                            variable = textscan(text,'%q','delimiter','=');
                            [Value] = textscan(variable{1}{2}(2:end-1),'%n%s');
                            results(aa).variable(bb).name = variable{1}{1};
                            results(aa).variable(bb).value = Value{1}*unit_scailer(Value{2});
                        end
                    end
                end
            end
            
            input_Line = fgetl(FID);
            
            count=0;
            while input_Line ~= -1;
                Values = textscan(input_Line,'%q','delimiter',',');
                count = count+1;
                for aa = 1:length(Values{1})
                    if isempty(Values{1}{aa})
                        results(aa).value(count)=nan;
                    else
                        results(aa).value(count)=str2num(Values{1}{aa});
                    end
                end
                input_Line = fgetl(FID);
                if rem(count,100) == 0
                    disp([num2str(count), ' lines compleated.']);
                end
            end
            fclose(FID);
            
            parameters(1).variable = results(2).variable;
            parameters(1).output(1).name = results(1).name;
            parameters(1).output(1).value = results(1).value;
            
            for aa = 2:length(results)
                same_variables_pnt = 1;
                while ~same(parameters(same_variables_pnt).variable,results(aa).variable);
                    same_variables_pnt = same_variables_pnt + 1;
                    if same_variables_pnt > length(parameters)
                        break
                    end
                end
                parameters(same_variables_pnt).variable = results(aa).variable;
                parameters(same_variables_pnt).output(1).name = results(1).name;
                parameters(same_variables_pnt).output(1).value = results(1).value*results(1).units;
                output_pnt = length(parameters(same_variables_pnt).output)+1;
                parameters(same_variables_pnt).output(output_pnt).name = results(aa).name;
                parameters(same_variables_pnt).output(output_pnt).value = results(aa).value*results(aa).units;
            end
            
            function [result] = same(first, second)
                
                result = 1;
                if length(first) == length(second)
                    for cc = 1:length(first)
                        if first(cc).value ~= second(cc).value
                            result = 0;
                        end
                    end
                else
                    result = 0;
                end
            end
            
            function [scailing] = unit_scailer(input)
                
                if isempty(input)
                    scailing = 1;
                elseif strcmp(input, 'cm')
                    scailing = 1e-2;
                elseif strcmp(input, 'mm')
                    scailing = 1e-3;
                elseif strcmp(input, 'um')
                    scailing = 1e-6;
                elseif strcmp(input, 'u')
                    scailing = 1e-6;
                elseif strcmp(input, 'mil')
                    scailing = 0;
                elseif strcmp(input, 'GHz')
                    scailing = 1e9;
                elseif strcmp(input, 'MHz')
                    scailing = 1e6;
                elseif strcmp(input, 'V')
                    scailing = 1;
                elseif strcmp(input, 'mV')
                    scailing = 1e-3;
                elseif strcmp(input, 'deg')
                    scailing = HFSS_Tools.degrees;
                elseif strcmp(input, 'rad')
                    scailing = 1;
                elseif strcmp(input, 'V_per_meter')
                    scailing = 1;
                elseif strcmp(input, 'eg') %To acount for matlab bug!!!
                    scailing = HFSS_Tools.degrees;
                else
                    scailing = 1;
                end
            end
        end
        
        
        function [ dB_voltage ] = dBv( voltage )
            %DBV dB from voltage
            
            dB_voltage = 20*log10(abs(voltage));
        end
        
        
        function [ dB_power ] = dBp( power )
            %DBP dB from power
            
            dB_power = 10*log10(abs(power));
        end
        
        
        function [ voltage ] = idBv( dB_voltage )
            %IDBV voltage magnitude from dB
            
            voltage = 10.^(dB_voltage/20);
        end
        
        
        function [ power ] = idBp( dB_power )
            %IDBV power from dB
            
            power = 10.^(dB_power./10);
        end
        
        
        function [ rad ] = degrees()
            %DEGREES multiplier to convert degrees to rad
            
            rad = pi/180;
        end
        
        
        function Free_Space_Velocity_Of_Light=c0()
            %c0 returns the free space velocity of light.
            
            Free_Space_Velocity_Of_Light=2.99792458e8;
        end
        
        
        function Permeability_of_Free_Space=u0()
            %u0 returns the permeability of free space.
            
            Permeability_of_Free_Space = 4*pi*1e-7;
        end
        
        
        function Permittivity_of_Free_Space=e0()
            %e0 returns the permittivity of free space.
            
            Permittivity_of_Free_Space = 1/(HFSS_Tools.u0()*HFSS_Tools.c0^2);
        end
        
        
        function Impedance_of_Free_Space=wave_impedance()
            %e0 returns the permittivity of free space.
            
            Impedance_of_Free_Space = HFSS_Tools.u0()*HFSS_Tools.c0;
        end
        
        
        function [ beam_pattern ] = rotate_axis(beam_pattern, roll, pitch, yaw)
            %GET_GAIN_THETA Find the realised theta gain of a beam pattern
            %  This function finds returns the theta polarised gain in the
            %  direction defined by theta and phi.
            
            for aa = 1:length(beam_pattern)
                [ beam_pattern_vector(aa) ] = HFSS_Tools.polar_to_vector(beam_pattern(aa));
                for bb = 1:numel(beam_pattern_vector(aa).direction_vector)
                    x = beam_pattern_vector(aa).direction_vector(bb).x;
                    y = beam_pattern_vector(aa).direction_vector(bb).y;
                    z = beam_pattern_vector(aa).direction_vector(bb).z;
                    [ x,y,z, ~ ] = Vector_Earth2Body( roll, pitch, yaw, x,y,z  );
                    beam_pattern_vector(aa).direction_vector(bb).x = x;
                    beam_pattern_vector(aa).direction_vector(bb).y = y;
                    beam_pattern_vector(aa).direction_vector(bb).z = z;
                    
                    x = beam_pattern_vector(aa).theta_pol_vector(bb).x;
                    y = beam_pattern_vector(aa).theta_pol_vector(bb).y;
                    z = beam_pattern_vector(aa).theta_pol_vector(bb).z;
                    [ x,y,z, ~ ] = Vector_Earth2Body( roll, pitch, yaw, x,y,z  );
                    beam_pattern_vector(aa).theta_pol_vector(bb).x = x;
                    beam_pattern_vector(aa).theta_pol_vector(bb).y = y;
                    beam_pattern_vector(aa).theta_pol_vector(bb).z = z;
                    
                    x = beam_pattern_vector(aa).phi_pol_vector(bb).x;
                    y = beam_pattern_vector(aa).phi_pol_vector(bb).y;
                    z = beam_pattern_vector(aa).phi_pol_vector(bb).z;
                    [ x,y,z, ~ ] = Vector_Earth2Body( roll, pitch, yaw, x,y,z  );
                    beam_pattern_vector(aa).phi_pol_vector(bb).x = x;
                    beam_pattern_vector(aa).phi_pol_vector(bb).y = y;
                    beam_pattern_vector(aa).phi_pol_vector(bb).z = z;
                end
                
                [ beam_pattern_polar(aa) ] = HFSS_Tools.vector_to_polar(beam_pattern_vector(aa), beam_pattern(aa).phi(:,1));
                x = beam_pattern_polar(aa).theta;
                y = beam_pattern_polar(aa).phi;
                xq = beam_pattern(aa).theta;
                yq = beam_pattern(aa).phi;
                
                v = beam_pattern_polar(aa).rE_phi;
                F = scatteredInterpolant(x',y',v','linear','nearest');
                beam_pattern(aa).realized_gain_phi = F(xq,yq);
                
                v = beam_pattern_polar(aa).rE_theta;
                F = scatteredInterpolant(x',y',v','linear','nearest');
                beam_pattern(aa).realized_gain_phi = F(xq,yq);
                
                v = beam_pattern_polar(aa).realized_gain_phi;
                F = scatteredInterpolant(x',y',v','linear','nearest');
                beam_pattern(aa).realized_gain_phi = F(xq,yq);
                
                v = beam_pattern_polar(aa).realized_gain_theta;
                F = scatteredInterpolant(x',y',v','linear','nearest');
                beam_pattern(aa).realized_gain_theta = F(xq,yq);
                
                v = beam_pattern_polar(aa).realized_gain_total;
                F = scatteredInterpolant(x',y',v','linear','nearest');
                beam_pattern(aa).realized_gain_total = F(xq,yq);
            end
        end
        
        
        function [ beam_pattern_polar ] = vector_to_polar(beam_pattern_vector, phi)
            %GET_GAIN_THETA Find the realised theta gain of a beam pattern
            %  This function finds returns the theta polarised gain in the
            %  direction defined by theta and phi.
            
            wrap_angle = pi/4;
            
            count = 0;
            for aa = 1:length(phi)
                count = count + 1;
                phi_out(count) = phi(aa);
                if abs(phi(aa)) > wrap_angle
                    count = count + 1;
                    if (phi(aa) > 0)
                        phi_out(count) = phi(aa)-pi;
                    else
                        phi_out(count) = phi(aa)+pi;
                    end
                end
            end
            
            count = 0;
            for bb = 1:numel(beam_pattern_vector.direction_vector)
                x = beam_pattern_vector.direction_vector(bb).x;
                y = beam_pattern_vector.direction_vector(bb).y;
                z = beam_pattern_vector.direction_vector(bb).z;
                [azimuth, elevation,~] = cart2sph(x,y,z);
                [ theta, phi ] = HFSS_Tools.az_el_to_theta_phi( azimuth, elevation );
                count = count + 1;
                beam_pattern_polar.theta(count) = theta;
                beam_pattern_polar.phi(count) = phi;
                
                
                [x,y,z] = sph2cart(azimuth,elevation+pi/2,1);
                theta_pol_vector_new = [x,y,z];
                
                [x,y,z] = sph2cart(azimuth+pi/2,0,1);
                phi_pol_vector_new = [x,y,z];
                
                x = beam_pattern_vector.theta_pol_vector(bb).x;
                y = beam_pattern_vector.theta_pol_vector(bb).y;
                z = beam_pattern_vector.theta_pol_vector(bb).z;
                theta_pol_vector_old = [x,y,z];
                
                x = beam_pattern_vector.phi_pol_vector(bb).x;
                y = beam_pattern_vector.phi_pol_vector(bb).y;
                z = beam_pattern_vector.phi_pol_vector(bb).z;
                phi_pol_vector_old = [x,y,z];
                
                angle_tt_pp = acos(dot(theta_pol_vector_new, theta_pol_vector_old));
                angle_tp = acos(dot(theta_pol_vector_new, phi_pol_vector_old));
                angle_pt = acos(dot(phi_pol_vector_new, theta_pol_vector_old));
                
                if angle_tp > angle_pt
                    pol_rotation_angle = angle_tt_pp;
                    rotation_matrix = [cos(pol_rotation_angle), -sin(pol_rotation_angle); sin(pol_rotation_angle), cos(pol_rotation_angle)];
                    Jones_vector = rotation_matrix * [beam_pattern_vector.rE_theta(bb); beam_pattern_vector.rE_phi(bb)];
                    rE_theta = Jones_vector(1);
                    rE_phi = Jones_vector(2);
                elseif angle_tp < angle_pt
                    pol_rotation_angle = -angle_tt_pp;
                    rotation_matrix = [cos(pol_rotation_angle), -sin(pol_rotation_angle); sin(pol_rotation_angle), cos(pol_rotation_angle)];
                    Jones_vector = rotation_matrix * [beam_pattern_vector.rE_theta(bb); beam_pattern_vector.rE_phi(bb)];
                    rE_theta = Jones_vector(1);
                    rE_phi = Jones_vector(2);
                else
                    rE_phi = beam_pattern_vector.rE_theta(bb);
                    rE_theta = beam_pattern_vector.rE_phi(bb);
                end
                beam_pattern_polar.rE_theta(count) = rE_theta;
                beam_pattern_polar.rE_phi(count) = rE_phi;
                beam_pattern_polar.realized_gain_theta(count) = abs(rE_theta.^2);
                beam_pattern_polar.realized_gain_phi(count) = abs(rE_phi.^2);
                beam_pattern_polar.realized_gain_total(count) = beam_pattern_vector.realized_gain_total(bb);
                
                if (abs(theta) < pi/720) || (abs(theta+pi) < pi/720) || (abs(theta-pi) < pi/720)
                    for cc = 1:length(phi_out)
                        count = count + 1;
                        beam_pattern_polar.theta(count) = theta;
                        beam_pattern_polar.phi(count) = phi_out(cc);
                        
                        
                        [x,y,z] = sph2cart(azimuth,elevation+pi/2,1);
                        theta_pol_vector_new = [x,y,z];
                        
                        [x,y,z] = sph2cart(azimuth+pi/2,0,1);
                        phi_pol_vector_new = [x,y,z];
                        
                        x = beam_pattern_vector.theta_pol_vector(bb).x;
                        y = beam_pattern_vector.theta_pol_vector(bb).y;
                        z = beam_pattern_vector.theta_pol_vector(bb).z;
                        theta_pol_vector_old = [x,y,z];
                        
                        x = beam_pattern_vector.phi_pol_vector(bb).x;
                        y = beam_pattern_vector.phi_pol_vector(bb).y;
                        z = beam_pattern_vector.phi_pol_vector(bb).z;
                        phi_pol_vector_old = [x,y,z];
                        
                        angle_tt_pp = acos(dot(theta_pol_vector_new, theta_pol_vector_old));
                        angle_tp = acos(dot(theta_pol_vector_new, phi_pol_vector_old));
                        angle_pt = acos(dot(phi_pol_vector_new, theta_pol_vector_old));
                        
                        if angle_tp > angle_pt
                            pol_rotation_angle = angle_tt_pp;
                            rotation_matrix = [cos(pol_rotation_angle), -sin(pol_rotation_angle); sin(pol_rotation_angle), cos(pol_rotation_angle)];
                            Jones_vector = rotation_matrix * [beam_pattern_vector.rE_theta(bb); beam_pattern_vector.rE_phi(bb)];
                            rE_theta = Jones_vector(1);
                            rE_phi = Jones_vector(2);
                        elseif angle_tp < angle_pt
                            pol_rotation_angle = -angle_tt_pp;
                            rotation_matrix = [cos(pol_rotation_angle), -sin(pol_rotation_angle); sin(pol_rotation_angle), cos(pol_rotation_angle)];
                            Jones_vector = rotation_matrix * [beam_pattern_vector.rE_theta(bb); beam_pattern_vector.rE_phi(bb)];
                            rE_theta = Jones_vector(1);
                            rE_phi = Jones_vector(2);
                        else
                            rE_phi = beam_pattern_vector.rE_theta(bb);
                            rE_theta = beam_pattern_vector.rE_phi(bb);
                        end
                        beam_pattern_polar.rE_theta(count) = rE_theta;
                        beam_pattern_polar.rE_phi(count) = rE_phi;
                        beam_pattern_polar.realized_gain_theta(count) = abs(rE_theta.^2);
                        beam_pattern_polar.realized_gain_phi(count) = abs(rE_phi.^2);
                        beam_pattern_polar.realized_gain_total(count) = beam_pattern_vector.realized_gain_total(bb);
                    end
                end
                
                if abs(phi) > wrap_angle
                    count = count + 1;
                    beam_pattern_polar.theta(count) = -theta;
                    if (phi > 0)
                        beam_pattern_polar.phi(count) = phi-pi;
                    else
                        beam_pattern_polar.phi(count) = phi+pi;
                    end
                    
                    
                    [x,y,z] = sph2cart(azimuth,elevation+pi/2,1);
                    theta_pol_vector_new = [x,y,z];
                    
                    [x,y,z] = sph2cart(azimuth+pi/2,0,1);
                    phi_pol_vector_new = [x,y,z];
                    
                    x = beam_pattern_vector.theta_pol_vector(bb).x;
                    y = beam_pattern_vector.theta_pol_vector(bb).y;
                    z = beam_pattern_vector.theta_pol_vector(bb).z;
                    theta_pol_vector_old = [x,y,z];
                    
                    x = beam_pattern_vector.phi_pol_vector(bb).x;
                    y = beam_pattern_vector.phi_pol_vector(bb).y;
                    z = beam_pattern_vector.phi_pol_vector(bb).z;
                    phi_pol_vector_old = [x,y,z];
                    
                    angle_tt_pp = acos(dot(theta_pol_vector_new, theta_pol_vector_old));
                    angle_tp = acos(dot(theta_pol_vector_new, phi_pol_vector_old));
                    angle_pt = acos(dot(phi_pol_vector_new, theta_pol_vector_old));
                    
                    if angle_tp > angle_pt
                        pol_rotation_angle = angle_tt_pp;
                        rotation_matrix = [cos(pol_rotation_angle), -sin(pol_rotation_angle); sin(pol_rotation_angle), cos(pol_rotation_angle)];
                        Jones_vector = rotation_matrix * [beam_pattern_vector.rE_theta(bb); beam_pattern_vector.rE_phi(bb)];
                        rE_theta = Jones_vector(1);
                        rE_phi = Jones_vector(2);
                    elseif angle_tp < angle_pt
                        pol_rotation_angle = -angle_tt_pp;
                        rotation_matrix = [cos(pol_rotation_angle), -sin(pol_rotation_angle); sin(pol_rotation_angle), cos(pol_rotation_angle)];
                        Jones_vector = rotation_matrix * [beam_pattern_vector.rE_theta(bb); beam_pattern_vector.rE_phi(bb)];
                        rE_theta = Jones_vector(1);
                        rE_phi = Jones_vector(2);
                    else
                        rE_phi = beam_pattern_vector.rE_theta(bb);
                        rE_theta = beam_pattern_vector.rE_phi(bb);
                    end
                    beam_pattern_polar.rE_theta(count) = rE_theta;
                    beam_pattern_polar.rE_phi(count) = rE_phi;
                    beam_pattern_polar.realized_gain_theta(count) = abs(rE_theta.^2);
                    beam_pattern_polar.realized_gain_phi(count) = abs(rE_phi.^2);
                    beam_pattern_polar.realized_gain_total(count) = beam_pattern_vector.realized_gain_total(bb);
                end
            end
        end
        
        
        function [ beam_pattern_vector ] = polar_to_vector(beam_pattern)
            %GET_GAIN_THETA Find the realised theta gain of a beam pattern
            %  This function finds returns the theta polarised gain in the
            %  direction defined by theta and phi.
            
            theta = beam_pattern.theta;
            phi = beam_pattern.phi;
            theta_0 = 1;
            theta_180 = 1;
            count = 0;
            for bb = 1:numel(theta)
                if ( (abs(theta(bb)) > pi/720) && (abs(theta(bb)+pi) > pi/720) && (abs(theta(bb)-pi) > pi/720) && (abs(phi(bb)+pi/2) > pi/720) ) || (abs(theta(bb)) <= pi/720 && theta_0) || (( (abs(theta(bb)+pi) <= pi/720) || (abs(theta(bb)-pi) <= pi/720) ) && theta_180)
                    if abs(theta(bb)) < pi/720
                        theta_0 = 0;
                    elseif (abs(theta(bb)+pi) < pi/720) || (abs(theta(bb)-pi) < pi/720)
                        theta_180 = 0;
                    end
                    count = count + 1;
                    [ azimuth, elevation ] = HFSS_Tools.theta_phi_to_az_el( theta(bb), phi(bb) );
                    [x,y,z] = sph2cart(azimuth,elevation,1);
                    beam_pattern_vector.direction_vector(count).x = x;
                    beam_pattern_vector.direction_vector(count).y = y;
                    beam_pattern_vector.direction_vector(count).z = z;
                    
                    [x,y,z] = sph2cart(azimuth,elevation+pi/2,1);
                    beam_pattern_vector.theta_pol_vector(count).x = x;
                    beam_pattern_vector.theta_pol_vector(count).y = y;
                    beam_pattern_vector.theta_pol_vector(count).z = z;
                    
                    [x,y,z] = sph2cart(azimuth+pi/2,0,1);
                    beam_pattern_vector.phi_pol_vector(count).x = x;
                    beam_pattern_vector.phi_pol_vector(count).y = y;
                    beam_pattern_vector.phi_pol_vector(count).z = z;
                    
                    beam_pattern_vector.rE_phi(count) = beam_pattern.rE_phi(bb);
                    beam_pattern_vector.rE_theta(count) = beam_pattern.rE_theta(bb);
                    beam_pattern_vector.realized_gain_phi(count) = beam_pattern.realized_gain_phi(bb);
                    beam_pattern_vector.realized_gain_theta(count) = beam_pattern.realized_gain_theta(bb);
                    beam_pattern_vector.realized_gain_total(count) = beam_pattern.realized_gain_total(bb);
                end
            end
        end
        
        
        function [ theta, phi ] = az_el_to_theta_phi( az, el )
            %UNWRAP_ANGLES Summary of this function goes here
            %   Detailed explanation goes here
            
            az = rem(az, 2*pi);
            if az < -pi
                az = az + 2*pi;
            elseif az > pi
                az = az - 2*pi;
            end
            theta = pi/2 - el;
            phi = az;
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
        end
        
        
        function [ az, el ] = theta_phi_to_az_el( theta, phi )
            %UNWRAP_ANGLES Summary of this function goes here
            %   Detailed explanation goes here
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            if theta < 0
                theta = -theta;
                if phi > 0
                    phi = phi - pi;
                else
                    phi = phi + pi;
                end
            end
            el = pi/2 - theta;
            az = phi;
        end
        
        
        function [ theta, phi ] = unwrap_angles( theta, phi )
            %UNWRAP_ANGLES Summary of this function goes here
            %   Detailed explanation goes here
            
            theta = rem(theta, 2*pi);
            if theta < -pi
                theta = theta + 2*pi;
            elseif theta > pi
                theta = theta - 2*pi;
            end
            phi = rem(phi, 2*pi);
            if phi < -pi
                phi = phi + 2*pi;
            elseif phi > pi
                phi = phi - 2*pi;
            end
            if phi < -pi/2
                phi = phi + pi;
                theta = -theta;
            elseif phi > pi/2
                phi = phi - pi;
                theta = -theta;
            end
        end
        
        
        function [ Body_X, Body_Y, Body_Z, DCM_EB ] = Vector_Earth2Body( Roll_A, Pitch_A, Yaw_A,  Earth_X, Earth_Y, Earth_Z  )
            %VECTOR_BODY2EARTH Summary of this function goes here
            %   Detailed explanation goes here
            
            
            DCM_EB = [cos(Pitch_A) * cos(Yaw_A), ...
                (cos(Pitch_A) * sin(Yaw_A)), ...
                - sin(Pitch_A); ...
                ...
                (-cos(Roll_A) * sin(Yaw_A) + sin(Roll_A) * sin(Pitch_A) * cos(Yaw_A)), ...
                (cos(Roll_A) * cos(Yaw_A) + sin(Roll_A) * sin(Pitch_A) * sin(Yaw_A)), ...
                sin(Roll_A) * cos(Pitch_A); ...
                ...
                (sin(Roll_A) * sin(Yaw_A) + cos(Roll_A) * sin(Pitch_A) * cos(Yaw_A)), ...
                (-sin(Roll_A) * cos(Yaw_A) + cos(Roll_A) * sin(Pitch_A) * sin(Yaw_A)), ...
                cos(Roll_A) * cos(Pitch_A)];
            
            Body_XYZ = DCM_EB * [Earth_X; Earth_Y; Earth_Z];
            
            Body_X = Body_XYZ(1);
            Body_Y = Body_XYZ(2);
            Body_Z = Body_XYZ(3);
        end
        
        
        function [ Earth_X, Earth_Y, Earth_Z, DCM_BE ] = Vector_Body2Earth( Roll_A, Pitch_A, Yaw_A, Body_X, Body_Y, Body_Z )
            %VECTOR_BODY2EARTH Summary of this function goes here
            %   Detailed explanation goes here
            
            % Earth_X = Body_X * cos(Pitch_A) * cos(Yaw_A) ...
            %     + Body_Y * (sin(Roll_A) * sin(Pitch_A) * cos(Yaw_A) - cos(Roll_A) * sin(Yaw_A)) ...
            %     + Body_Z * (cos(Roll_A) * sin(Pitch_A) * cos(Yaw_A) + sin(Roll_A) * sin(Yaw_A));
            %
            % Earth_Y = Body_X * cos(Pitch_A) * sin(Yaw_A) ...
            %     + Body_Y * (sin(Roll_A) * sin(Pitch_A) * sin(Yaw_A) + cos(Roll_A) * cos(Yaw_A)) ...
            %     + Body_Z * (cos(Roll_A) * sin(Pitch_A) * sin(Yaw_A) - sin(Roll_A) * cos(Yaw_A));
            %
            % Earth_Z = - Body_X * sin(Pitch_A) ...
            %     + Body_Y * sin(Roll_A) * cos(Pitch_A) ...
            %     + Body_Z * cos(Roll_A) * cos(Pitch_A);
            
            DCM_BE = [cos(Pitch_A) * cos(Yaw_A) ...
                (sin(Roll_A) * sin(Pitch_A) * cos(Yaw_A) - cos(Roll_A) * sin(Yaw_A)) ...
                (cos(Roll_A) * sin(Pitch_A) * cos(Yaw_A) + sin(Roll_A) * sin(Yaw_A)); ...
                ...
                cos(Pitch_A) * sin(Yaw_A) ...
                (sin(Roll_A) * sin(Pitch_A) * sin(Yaw_A) + cos(Roll_A) * cos(Yaw_A)) ...
                (cos(Roll_A) * sin(Pitch_A) * sin(Yaw_A) - sin(Roll_A) * cos(Yaw_A)); ...
                ...
                - sin(Pitch_A) ...
                sin(Roll_A) * cos(Pitch_A) ...
                cos(Roll_A) * cos(Pitch_A)];
            
            Earth_XYZ = DCM_BE * [Body_X; Body_Y; Body_Z];
            
            Earth_X = Earth_XYZ(1);
            Earth_Y = Earth_XYZ(2);
            Earth_Z = Earth_XYZ(3);
        end
        
        
        function [ realized_gain_theta, frequency ] = get_gain_theta(beam_pattern, theta, phi)
            %GET_GAIN_THETA Find the realised theta gain of a beam pattern
            %  This function finds returns the theta polarised gain in the
            %  direction defined by theta and phi.
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            
            realized_gain_theta = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                X = beam_pattern(aa).theta;
                Y = beam_pattern(aa).phi;
                Z = beam_pattern(aa).realized_gain_theta;
                realized_gain_theta(aa) = interp2(X,Y,Z,theta,phi);
                frequency(aa) = beam_pattern(aa).frequency;
            end
        end
        
        
        function [ realised_gain_phi, frequency ] = get_gain_phi(beam_pattern, theta, phi)
            %GET_GAIN_PHI Find the realised phi gain of a beam pattern
            %  This function finds returns the phi polarised gain in the
            %  direction defined by theta and phi.
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            
            realised_gain_phi = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                X = beam_pattern(aa).theta;
                Y = beam_pattern(aa).phi;
                Z = beam_pattern(aa).realized_gain_phi;
                realised_gain_phi(aa) = interp2(X,Y,Z,theta,phi);
                frequency(aa) = beam_pattern(aa).frequency;
            end
        end
        
        
        function [ realized_gain_total, frequency ] = get_gain_total(beam_pattern, theta, phi)
            %GET_GAIN_TOTAL Find the realised gain of a beam pattern
            %  This function finds returns the unpolarised gain in the
            %  direction defined by theta and phi.
            
            [ theta, phi ] = HFSS_Tools.unwrap_angles( theta, phi );
            
            realized_gain_total = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                X = beam_pattern(aa).theta;
                Y = beam_pattern(aa).phi;
                Z = beam_pattern(aa).realized_gain_total;
                realized_gain_total(aa) = interp2(X,Y,Z,theta,phi);
                frequency(aa) = beam_pattern(aa).frequency;
            end
        end
        
        
        function [ realized_gain_theta, theta, phi, frequency ] = get_max_gain_theta(beam_pattern)
            %GET_MAX_GAIN_THETA Find the maximum realised theta gain of a beam pattern
            %  This function finds the direction of the maximum realied
            %  gain polarised in the theta direction. It outputs the maximum
            %  gain and the angles that it found that gain, theta and phi.
            
            realized_gain_theta = ones(length(beam_pattern),1);
            theta = ones(length(beam_pattern),1);
            phi = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                [~, index] = max(beam_pattern(aa).realized_gain_theta(:));
                realized_gain_theta(aa) = beam_pattern(aa).realized_gain_theta(index);
                theta(aa) = beam_pattern(aa).theta(index);
                phi(aa) = beam_pattern(aa).phi(index);
                frequency(aa) = beam_pattern(aa).frequency;
            end
        end
        
        
        function [ realized_gain_phi, theta, phi, frequency ] = get_max_gain_phi(beam_pattern)
            %GET_MAX_GAIN_PHI Find the maximum realised phi gain of a beam pattern
            %  This function finds the direction of the maximum realied
            %  gain polarised in the phi direction. It outputs the maximum
            %  gain and the angles that it found that gain, theta and phi.
            
            realized_gain_phi = ones(length(beam_pattern),1);
            theta = ones(length(beam_pattern),1);
            phi = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                [~, index] = max(beam_pattern(aa).realized_gain_phi(:));
                realized_gain_phi(aa) = beam_pattern(aa).realized_gain_phi(index);
                theta(aa) = beam_pattern(aa).theta(index);
                phi(aa) = beam_pattern(aa).phi(index);
                frequency(aa) = beam_pattern(aa).frequency;
            end
        end
        
        
        function [ realised_gain, theta, phi, frequency ] = get_max_gain_total(beam_pattern)
            %GET_MAX_GAIN_TOTAL Find the maximum realised total gain of a beam pattern
            %  This function finds the direction of the maximum realied
            %  gain. It outputs the maximum gain and the angles that it
            %  found that gain, theta and phi.
            
            realised_gain = ones(length(beam_pattern),1);
            theta = ones(length(beam_pattern),1);
            phi = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                [~, index] = max(beam_pattern(aa).realized_gain_total(:));
                realised_gain(aa) = beam_pattern(aa).realized_gain_total(index);
                theta(aa) = beam_pattern(aa).theta(index);
                phi(aa) = beam_pattern(aa).phi(index);
                frequency(aa) = beam_pattern(aa).frequency;
            end
        end
        
        
        function [ theta_beam_width, phi_beam_width, frequency ] = get_beam_width(beam_pattern)
            %GET_MAX_GAIN_THETA Find the maximum realised theta gain of a beam pattern
            %  This function finds the direction of the maximum realied
            %  gain polarised in the theta direction. It outputs the maximum
            %  gain and the angles that it found that gain, theta and phi.
            
            realised_gain = ones(length(beam_pattern),1);
            theta = ones(length(beam_pattern),1);
            phi = ones(length(beam_pattern),1);
            theta_beam_width = zeros(length(beam_pattern),1);
            phi_beam_width = zeros(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                [~, index] = max(beam_pattern(aa).realized_gain_total(:));
                realised_gain(aa) = beam_pattern(aa).realized_gain_total(index);
                theta(aa) = beam_pattern(aa).theta(index);
                phi(aa) = beam_pattern(aa).phi(index);
                frequency(aa) = beam_pattern(aa).frequency;
                
                [ az, el ] = HFSS_Tools.theta_phi_to_az_el( theta(aa), phi(aa) );
                [ beam_pattern_temp ] = HFSS_Tools.rotate_axis(beam_pattern(aa), 0, -el, az);
                
                level = 0.5*realised_gain(aa);
                contour_temp = contourc(beam_pattern_temp.theta(1,:), beam_pattern_temp.phi(:,1).', beam_pattern_temp.realized_gain_total, [level, level]);
                contour_pnt = 0;
                pnt = 1;
                while pnt < length(contour_temp)
                    contour_pnt = contour_pnt + 1;
                    contour_3db(contour_pnt).data = contour_temp(:,pnt+1:contour_temp(2,pnt)+pnt);
                    pnt = contour_temp(2,pnt)+pnt+1;
                end
                theta_beam_width(aa) = 0;
                phi_beam_width(aa) = 0;
                for bb = 1:length(contour_3db)
                    theta_min = min(contour_3db(bb).data(1,:));
                    theta_max = max(contour_3db(bb).data(1,:));
                    phi_min = min(contour_3db(bb).data(2,:));
                    phi_max = max(contour_3db(bb).data(2,:));
                    if theta_min < pi/2 && theta_max > pi/2 && ...
                            phi_min < 0 && phi_max > 0
                        theta_beam_width(aa) = theta_max - theta_min;
                        phi_beam_width(aa) = phi_max - phi_min;
                    end
                end
            end
        end
        
        
        function [ side_lobe_gain, theta, phi, frequency ] = get_side_lobe(beam_pattern)
            %GET_MAX_GAIN_THETA Find the maximum realised theta gain of a beam pattern
            %  This function finds the direction of the maximum realied
            %  gain polarised in the theta direction. It outputs the maximum
            %  gain and the angles that it found that gain, theta and phi.
            
            %mask = poly2mask(x,y,12,12)
            %inpolygon
            
            side_lobe_gain = zeros(length(beam_pattern),1);
            theta = ones(length(beam_pattern),1);
            phi = ones(length(beam_pattern),1);
            frequency = ones(length(beam_pattern),1);
            for aa = 1:length(beam_pattern)
                [size_x, size_y] = size(beam_pattern(aa).realized_gain_total);
                [mesh_x,mesh_y] = meshgrid(1:size_y,1:size_x);
                [~, index] = max(beam_pattern(aa).realized_gain_total(:));
                [pnt_theta,pnt_phi] = ind2sub([size_x, size_y],index);
                realised_gain = beam_pattern(aa).realized_gain_total(index);
                theta(aa) = beam_pattern(aa).theta(index);
                phi(aa) = beam_pattern(aa).phi(index);
                frequency(aa) = beam_pattern(aa).frequency;
                
                [ az, el ] = HFSS_Tools.theta_phi_to_az_el( theta(aa), phi(aa) );
                [ beam_pattern_temp ] = HFSS_Tools.rotate_axis(beam_pattern(aa), 0, -el, -az);
                
                [~, index] = max(beam_pattern_temp.realized_gain_total(:));
                [pnt_theta,pnt_phi] = ind2sub([size_x, size_y],index);
                
                found_sidelobe = 0;
                level = realised_gain;
                while ~found_sidelobe
                    level = 0.5*level;
                    contour_temp = contourc(mesh_x(1,:), mesh_y(:,1).', beam_pattern_temp.realized_gain_total, [level, level]);
                    contour_pnt = 0;
                    pnt = 1;
                    contour_3db = [];
                    while pnt < length(contour_temp)
                        contour_pnt = contour_pnt + 1;
                        contour_3db(contour_pnt).data = contour_temp(:,pnt+1:contour_temp(2,pnt)+pnt);
                        pnt = contour_temp(2,pnt)+pnt+1;
                    end
                    if length(contour_3db) > 1
                        for bb = 1:length(contour_3db)
                            data = contour_3db(bb).data;
                            if abs(sum(wrapToPi(diff(angle((data(1,:) - pnt_phi) + (data(2,:) - pnt_theta)*1j))))) > pi;
%                                 figure
%                                 surf(beam_pattern(aa).phi/HFSS_Tools.degrees, beam_pattern(aa).theta/HFSS_Tools.degrees, dBp(beam_pattern_temp.realized_gain_total),'EdgeColor','none');
%                                 shading interp;
%                                 figure
%                                 hold on;
%                                 contour(beam_pattern(aa).phi/HFSS_Tools.degrees, beam_pattern(aa).theta/HFSS_Tools.degrees, dBp(beam_pattern_temp.realized_gain_total), dBp([level, level]));
%                                 contour3(beam_pattern(aa).phi/HFSS_Tools.degrees, beam_pattern(aa).theta/HFSS_Tools.degrees, dBp(beam_pattern_temp.realized_gain_total)-20, 50);
                                data_round = round(data);
                                mask = poly2mask(data_round(1,:),data_round(2,:),size_x,size_y);
                                realized_gain_total = beam_pattern_temp.realized_gain_total.*(1-mask);
                                [~, index] = max(realized_gain_total(:));
                                if side_lobe_gain(aa) < beam_pattern(aa).realized_gain_total(index)
                                    side_lobe_gain(aa) = beam_pattern(aa).realized_gain_total(index);
                                    theta(aa) = beam_pattern(aa).theta(index);
                                    phi(aa) = beam_pattern(aa).phi(index);
                                    frequency(aa) = beam_pattern(aa).frequency;
                                end
                                found_sidelobe = 1;
                            end
                        end
                    end
                end
            end
        end
        
        
        function [ pointer ] = plot_spherical_pattern( phi_grid, theta_grid, plot_mag, min_plot_mag )
            %PRESENT_ANT_SPH Summary of this function goes here
            %   Detailed explanation goes here
            
            switch nargin
                case 3
                    min_plot_mag = min(min(plot_mag));
                case 4
                    plot_mag = plot_mag - min_plot_mag;
                    plot_mag = plot_mag.*(plot_mag>0)+min_plot_mag;
                otherwise
                    disp('Wrong number of inputs');
            end
            [x,y,z] = sph2cart(phi_grid,pi/2-theta_grid,plot_mag-min_plot_mag);
            pointer = surf(x,y,z,plot_mag);
            MaxX = max(max(x));
            MaxY = max(max(y));
            MaxZ = max(max(z));
            line([0,1.1*MaxX],[0,0],[0,0],'Color','k');
            text(1.2*MaxX,0,0,'X')
            line([0,0],[0,1.1*MaxY],[0,0],'Color','k');
            text(0,1.2*MaxY,0,'Y')
            line([0,0],[0,0],[0,1.1*MaxZ],'Color','k');
            text(0,0,1.2*MaxZ,'Z')
            
            Line_Angle = 0:1*HFSS_Tools.degrees:10*HFSS_Tools.degrees;
            line(cos(Line_Angle)*1.1*MaxX,sin(Line_Angle)*1.1*MaxX,zeros(1,length(Line_Angle)),'Color','k');
            text(cos(12*HFSS_Tools.degrees)*1.1*MaxX,sin(12*HFSS_Tools.degrees)*1.1*MaxX,0,'\phi')
            line(sin(Line_Angle)*1.1*MaxZ,zeros(1,length(Line_Angle)),cos(Line_Angle)*1.1*MaxZ,'Color','k');
            text(sin(12*HFSS_Tools.degrees)*1.1*MaxZ,0,cos(12*HFSS_Tools.degrees)*1.1*MaxZ,'\theta')
            colorbar;
            axis equal;
            shading interp;
            axis off;
            colormap('jet');
        end
        
        
        function [ pointer ] = plot_surf_pattern( phi_grid, theta_grid, plot_mag, min_plot_mag )
            %PRESENT_ANT_SPH Summary of this function goes here
            %   Detailed explanation goes here
            
            pointer = surf(phi_grid/HFSS_Tools.degrees, theta_grid/HFSS_Tools.degrees, plot_mag,'EdgeColor','none');
            colormap('jet');
            shading interp;
            
            switch nargin
                case 3
                case 4
                    caxis([min_plot_mag,max(max(plot_mag))]);
                otherwise
                    disp('Wrong number of inputs');
            end
            colorbar;
            xlabel('phi (deg)');
            ylabel('theta (deg)');
        end
        
        
        function MrSmith(draw_admitance)
            %MRSMITH plot the smith chart
            %   This function plots a smith chart. To plot the impedance
            %   smith chart set draw_admitance = 0, to plot the admitance smith
            %   chart set draw_admitance = 1, default is 0.
            
            % Default values
            coarsfactor = 100;
            rvalues = [0 0.2 0.5 1 2 5];
            xvalues = [0.2 0.5 1 2 5];
            
            if (nargin == 0)
                draw_admitance = 0;
            end
            
            %We dont want any axis and equal x and y ratio
            axis('equal')
            axis('off')
            hold on;
            % Draw the horisontal line
            plot([-1 1], [0 0], 'k');
            % Draw the r-circles
            for r = rvalues
                theta=linspace(-pi, pi, coarsfactor);
                z = (r/(r+1) + 1/(r+1)*exp(1j*theta));
                if( draw_admitance )
                    z = -z;
                end
                plot(z, 'k');
                
                xpos = ZedToGamma(r);
                if(draw_admitance)
                    xpos = -xpos;
                end;
                handle=text(xpos, 0, num2str(r),'FontSize', get(gca,'FontSize'));
                set(handle, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
            end
            
            % Draw the x-circles
            tmp=linspace(0,sqrt(100), coarsfactor);
            r = tmp.^2;
            r = [0 r];
            for x = abs(xvalues)
                g = ZedToGamma(r+1j*x*ones(1, coarsfactor+1));
                if( draw_admitance )
                    g = -g;
                end
                plot(g, 'k');
                plot(conj(g), 'k');
                
                xpos = real(ZedToGamma(1j*x));
                ypos = imag(ZedToGamma(1j*x));
                if(draw_admitance)
                    xpos = -xpos;
                end;
                handle=text([xpos xpos], [ypos -ypos], [' j' num2str(x) ; '-j' num2str(x)],'FontSize', get(gca,'FontSize'));
                set(handle(1),'VerticalAlignment', 'bottom');
                set(handle(2),'VerticalAlignment', 'top');
                if xpos == 0
                    set(handle, 'Horizontalalignment', 'center');
                elseif xpos < 0
                    set(handle, 'Horizontalalignment', 'right');
                end
            end
            
            function [Gamma] = ZedToGamma (Zed)
                % Returns the Gamma(Reflection factor) values for the given Z (impedance)
                Gamma = (Zed - 1)./(Zed + 1);
            end
        end
        
        
    end
    
    
    methods(Access = private)
        %% HFSS Script Methods
        
        function [ hfss_script ] = HFSS_Open(obj, hfss_script )
            %HFSS_OPEN generates the VB script to select a project and design in HFSS
            
            hfss_script = [hfss_script, 'Dim oAnsoftApp \n'];
            hfss_script = [hfss_script, 'Dim oDesktop \n'];
            hfss_script = [hfss_script, 'Dim oProject \n'];
            hfss_script = [hfss_script, 'Dim oDesign \n'];
            hfss_script = [hfss_script, 'Dim oEditor \n'];
            hfss_script = [hfss_script, 'Dim oModule \n'];
            hfss_script = [hfss_script, 'Set oAnsoftApp = CreateObject("Ansoft.ElectronicsDesktop") \n'];
            hfss_script = [hfss_script, 'Set oDesktop = oAnsoftApp.GetAppDesktop() \n'];
            hfss_script = [hfss_script, 'oDesktop.RestoreWindow \n'];
            hfss_script = [hfss_script, 'Set oProject = oDesktop.SetActiveProject("', obj.project_name, '") \n'];
            hfss_script = [hfss_script, 'Set oDesign = oProject.SetActiveDesign("', obj.active_design, '") \n'];
        end
        
        
        function [ hfss_script ] = HFSS_Set_Mode_Sources(obj, hfss_script )
            %HFSS_SET_MODE_SOURCES generates the VB script to set sorces in HFSS
            
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("Solutions") \n'];
            hfss_script = [hfss_script, 'oModule.EditSources "TotalFields", Array("NAME:Names"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', obj.sources(aa).name, '"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Modes"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', int2str(obj.sources(aa).modes), '"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Magnitudes"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', int2str(abs(obj.sources(aa).excitation_E)), 'W"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Phases"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', int2str(angle(obj.sources(aa).excitation_E)), 'rad"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  false \n'];
        end
        
        
        function [ hfss_script ] = HFSS_Set_Terminal_Sources(obj, hfss_script )
            %HFSS_SET_TERMINAL_SOURCES generates the VB script to set sorces in HFSS
            
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("Solutions") \n'];
            hfss_script = [hfss_script, 'oModule.EditSources "TotalFields",  _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Names"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', obj.sources(aa).name, '"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Terminals"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', ', int2str(obj.sources(aa).modes)];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Magnitudes"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', int2str(abs(obj.sources(aa).excitation_E)), 'W"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Phases"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', "', int2str(angle(obj.sources(aa).excitation_E)), 'rad"'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Terminated"'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, ', false'];
            end
            hfss_script = [hfss_script, '), _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:Impedances"), false, true \n'];
        end
        
        
        function [ hfss_script ] = DESIGNER_Set_Sources(obj, hfss_script )
            %HFSS_SET_TERMINAL_SOURCES generates the VB script to set sorces in HFSS
            
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("Excitations") \n'];
            hfss_script = [hfss_script, 'oModule.EditExcitations Array("NAME:Excitations",  _ \n'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, '   Array('];
                hfss_script = [hfss_script, '"NAME:', obj.sources(aa).name, '"'];
                hfss_script = [hfss_script, ', "', int2str(angle(obj.sources(aa).excitation_E)), 'rad"'];
                hfss_script = [hfss_script, ', "', int2str(abs(obj.sources(aa).excitation_E)), 'V"),  _ \n'];
            end
            hfss_script = [hfss_script(1:end-7), '),  _ \n'];
            hfss_script = [hfss_script, '   Array("NAME:PostProcess",  _ \n'];
            for aa = 1:obj.source_number;
                hfss_script = [hfss_script, '   Array('];
                hfss_script = [hfss_script, '"NAME:', obj.sources(aa).name, '", true, "0mm", "50ohm + 0i ohm"),  _ \n'];
            end
            hfss_script = [hfss_script(1:end-7), '),  _ \n'];
            hfss_script = [hfss_script, '  "GapSource,IncidentVoltage",  _ \n'];
            hfss_script = [hfss_script, '  Array("NAME:PushExParamsBlock", "IsTransient:=", false, "StartTime:=", "0s",  _ \n'];
            hfss_script = [hfss_script, '  "StopTime:=", "10s", "MaxHarmonics:=", 100, "WinType:=", 0, "WidthPercentage:=", 100,  _ \n'];
            hfss_script = [hfss_script, '  "KaiserParam:=", 0) \n'];
        end
        
        
        function [ hfss_script ] = HFSS_Add_Rad_Sphere(obj, hfss_script )
            %HFSS_ADD_RAD_SPHERE generates the VB script to add a radiation sphere in HFSS
            
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("RadField") \n'];
            hfss_script = [hfss_script, 'oModule.InsertFarFieldSphereSetup Array("NAME:Matlab Sphere", "UseCustomRadiationSurface:=",  _ \n'];
            hfss_script = [hfss_script, '  false,  _ \n'];
            hfss_script = [hfss_script, '  "ThetaStart:=", "-180deg", "ThetaStop:=", "180deg", "ThetaStep:=", "', num2str(obj.theta_step), 'deg",  _ \n'];
            hfss_script = [hfss_script, '  "PhiStart:=", "-90deg", "PhiStop:=", "90deg", "PhiStep:=", "', num2str(obj.phi_step), 'deg",  _ \n'];
            
            if isempty(obj.fields_axis)
                hfss_script = [hfss_script, '  "UseLocalCS:=", false) \n'];
                hfss_script = [hfss_script, ' \n'];
            else
                hfss_script = [hfss_script, '  "UseLocalCS:=", true, "CoordSystem:=", "', obj.fields_axis, '") \n'];
                hfss_script = [hfss_script, ' \n'];
            end
        end
        
        
        function [ hfss_script ] = HFSS_Add_Rad_Report(obj, hfss_script )
            %HFSS_ADD_RAD_REPORT generates the VB script to add a radiation report in HFSS
            
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("ReportSetup") \n'];
            hfss_script = [hfss_script, 'oModule.CreateReport "Matlab Table", "Far Fields", "Data Table", "'];
            if isempty(obj.solution_sweep_fields)
                hfss_script = [hfss_script, obj.solution_setup];
                hfss_script = [hfss_script, ' : LastAdaptive",  _ \n'];
            else
                hfss_script = [hfss_script, obj.solution_setup];
                hfss_script = [hfss_script, ' : '];
                hfss_script = [hfss_script, obj.solution_sweep_fields];
                hfss_script = [hfss_script, '",  _ \n'];
            end
            hfss_script = [hfss_script, '  Array("Context:=", "Matlab Sphere"),  _ \n'];
            hfss_script = [hfss_script, '  Array("Theta:=", Array( "All"), "Phi:=", Array("All"), "Freq:=", Array("All")),  _ \n'];
            hfss_script = [hfss_script, '  Array("X Component:=", "Theta", "Y Component:=", Array( "rEPhi", "rETheta")), Array() \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
        
        function [ hfss_script ] = DESIGNER_Add_Rad_Report(obj, hfss_script )
            %HFSS_ADD_RAD_REPORT generates the VB script to add a radiation report in HFSS
            
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("ReportSetup") \n'];
            hfss_script = [hfss_script, 'oModule.CreateReport "Matlab Table", "Far Field", "Data Table",  _ \n'];
            if isempty(obj.solution_sweep_fields)
            else
                hfss_script = [hfss_script, '"', obj.solution_setup, ' : ', obj.solution_sweep_fields, '", Array("NAME:Context", "SimValueContext:=", Array(3, 0, 2,  _ \n'];
            end
            hfss_script = [hfss_script, '0, false, false, -1, 1, 0, 1, 1, "", 0, 0, "EnsDiffPairKey", false, "0", "IDIID",  _ \n'];
            hfss_script = [hfss_script, 'false, "1")), Array(  _ \n'];
            
            hfss_script = [hfss_script, '"Theta:=", Array("All"), "OverridingValues:=", Array(  _ \n'];
            for Theta_angle = -180:obj.theta_step:180-obj.theta_step
                hfss_script = [hfss_script, '"', num2str(Theta_angle), 'deg",'];
            end
            hfss_script = [hfss_script, ' "180deg"),  _ \n'];
            hfss_script = [hfss_script, '"Phi:=", Array("All"), "OverridingValues:=", Array(  _ \n'];
            for Phi_angle = -90:obj.phi_step:90-obj.phi_step
                hfss_script = [hfss_script, '"', num2str(Phi_angle), 'deg",'];
            end
            hfss_script = [hfss_script, ' "90deg"),  _ \n'];
            
            hfss_script = [hfss_script, '"F:=", Array("All")),  _ \n'];
            hfss_script = [hfss_script, 'Array("X Component:=", "Theta", "Y Component:=",  _ \n'];
            hfss_script = [hfss_script, 'Array( "Ephi", "Etheta")), Array() \n'];
        end
        
        
        function [ hfss_script ] = HFSS_S_Parameters(obj, hfss_script, file_name )
            %HFSS_S_PARAMETERS  generates the VB script to output s-parameters in HFSS
            
            file_name = strrep(file_name, '\', '/');
            
            hfss_script = [hfss_script, 'oDesign.ExportNetworkData '];
            hfss_script = [hfss_script, '"", _\n'];
            if isempty(obj.solution_sweep_freq)
                hfss_script = [hfss_script, sprintf(' Array("%s:LastAdaptive"), ', obj.solution_setup)];
            else
                hfss_script = [hfss_script, sprintf(' Array("%s:%s"), ', obj.solution_setup, obj.solution_sweep_freq)];
            end
            hfss_script = [hfss_script, sprintf('3, "%s", Array("All"), true, 50, "S", -1, 0, 15, true \n', file_name)];
        end
        
        
        function [ hfss_script ] = HFSS_Simulate_freq(obj, hfss_script )
            %HFSS_SIMULATE_FREQ generates the VB script to run the frequency simulation in HFSS
            
            if isempty(obj.solution_sweep_freq)
                hfss_script = [hfss_script, 'oDesign.Analyze "', obj.solution_setup, '" \n'];
            else
                hfss_script = [hfss_script, 'oDesign.Analyze "', obj.solution_setup, ' : ', obj.solution_sweep_freq, '" \n'];
            end
            hfss_script = [hfss_script, 'oProject.Save \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
        
        function [ hfss_script ] = HFSS_Simulate_fields(obj, hfss_script )
            %HFSS_SIMULATE_FIELDS generates the VB script to run the fields simulation in HFSS
            
            if isempty(obj.solution_sweep_fields)
                hfss_script = [hfss_script, 'oDesign.Analyze "', obj.solution_setup, '" \n'];
            else
                hfss_script = [hfss_script, 'oDesign.Analyze "', obj.solution_setup, ' : ', obj.solution_sweep_fields, '" \n'];
            end
            hfss_script = [hfss_script, 'oProject.Save \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
    end
    
    
    methods(Access = private, Static = true)
        %% HFSS Script Methods
        function [ hfss_script ] = HFSS_Save_Rad_Report(hfss_script, file_name )
            %HFSS_SAVE_RAD_REPORT generates the VB script to export the fields report in HFSS
            
            file_name = strrep(file_name, '\', '/');
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("ReportSetup") \n'];
            hfss_script = [hfss_script, 'oModule.ExportToFile "Matlab Table", "', file_name, '.csv" \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
        function [ hfss_script ] = HFSS_Remove_Matlab(hfss_script )
            %HFSS_REMOVE_MATLAB generates the VB script to remove the far field sphere and fields report in HFSS
            
            hfss_script = [hfss_script, 'oModule.DeleteReports Array("Matlab Table") \n'];
            hfss_script = [hfss_script, ' \n'];
            hfss_script = [hfss_script, 'Set oModule = oDesign.GetModule("RadField") \n'];
            hfss_script = [hfss_script, 'oModule.DeleteFarFieldSetup Array("Matlab Sphere") \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
        function [ hfss_script ] = DESIGNER_Remove_Matlab(hfss_script )
            %HFSS_REMOVE_MATLAB generates the VB script to remove the far field sphere and fields report in HFSS
            
            hfss_script = [hfss_script, 'oModule.DeleteReports Array("Matlab Table") \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
        function [ hfss_script ] = HFSS_Save_Parameters(hfss_script, param_name, param_value, param_units )
            %HFSS_SAVE_PARAMETERS generates the VB script to set parameters in HFSS
            
            if strcmp('none', param_units);
                param_units = '';
            end
            hfss_script = [hfss_script, 'oDesign.ChangeProperty Array("NAME:AllTabs", _\n'];
            hfss_script = [hfss_script, 'Array("NAME:LocalVariableTab", _\n'];
            hfss_script = [hfss_script, 'Array("NAME:PropServers", "LocalVariables"), _\n'];
            hfss_script = [hfss_script, 'Array("NAME:ChangedProps", _\n'];
            hfss_script = [hfss_script, 'Array("NAME:', param_name, '", "Value:=", "', num2str(param_value, 10), param_units, '")))) \n'];
            hfss_script = [hfss_script, ' \n'];
        end
        
    end
    
end


