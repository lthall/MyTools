clc
clear all
close all

% olympic = HFSS_Tools('TC1_ANT_Layout_090215_EXP', 'full_01_H_shortedV', 34);
% olympic = olympic.set_solution_setup('Setup1');
% olympic = olympic.set_solution_sweep_freq('Sweep');
% olympic = olympic.set_solution_sweep_fields('Sweep2');
% olympic = olympic.set_source_type('Modal');
% olympic = olympic.set_angle_resolution(5,5);
% 
% % we must manually match the order in HFSS.
% olympic = olympic.set_source_name( 1, '2v');
% olympic = olympic.set_source_name( 2, '2h');
% olympic = olympic.set_source_name( 3, '3v');
% olympic = olympic.set_source_name( 4, '3h');
% olympic = olympic.set_source_name( 5, '4v');
% olympic = olympic.set_source_name( 6, '4h');
% olympic = olympic.set_source_name( 7, '5v');
% olympic = olympic.set_source_name( 8, '5h');
% olympic = olympic.set_source_name( 9, '6v');
% olympic = olympic.set_source_name( 10, '6h');
% olympic = olympic.set_source_name( 11, '7v');
% olympic = olympic.set_source_name( 12, '7h');
% olympic = olympic.set_source_name( 13, '8v');
% olympic = olympic.set_source_name( 14, '8h');
% olympic = olympic.set_source_name( 15, '9h');
% olympic = olympic.set_source_name( 16, '10v');
% olympic = olympic.set_source_name( 17, '10h');
% olympic = olympic.set_source_name( 18, '11v');
% olympic = olympic.set_source_name( 19, '11h');
% olympic = olympic.set_source_name( 20, '12v');
% olympic = olympic.set_source_name( 21, '12h');
% olympic = olympic.set_source_name( 22, '13v');
% olympic = olympic.set_source_name( 23, '13h');
% olympic = olympic.set_source_name( 24, '14v');
% olympic = olympic.set_source_name( 25, '14h');
% olympic = olympic.set_source_name( 26, '15v');
% olympic = olympic.set_source_name( 27, '15h');
% olympic = olympic.set_source_name( 28, '16v');
% olympic = olympic.set_source_name( 29, '16h');
% olympic = olympic.set_source_name( 30, 'c_v');
% olympic = olympic.set_source_name( 31, 'c_h');
% olympic = olympic.set_source_name( 32, '9v');
% olympic = olympic.set_source_name( 33, '1v');
% olympic = olympic.set_source_name( 34, '1h');
% 
% % % olympic = olympic.simulate_freq;
% % % olympic = olympic.simulate_fields;
% 
% %% S-Parameters
% 
% %olympic = olympic.DESIGNER_export_S_parameters;
% %olympic = olympic.HFSS_export_S_parameters;
% olympic = olympic.import_S_parameters;
% 
% %% Fields
% 
% %olympic = olympic.HFSS_export_fields_data;
% olympic = olympic.import_fields_data;
% 
% %%
% 
% save('olympicDesigner', 'olympic');
% break;

%%

load('olympicDesigner');

%%

olympic = olympic.set_source_groups([21, 29, 27, 17, 12, 10, 6, 14, 8, 2], 1);
olympic = olympic.set_source_groups([15, 19, 34, 4], 2);
olympic = olympic.set_source_groups([25, 23], 3);

olympic.print_S_param(1);

% Plot S parameters

[frequency, S] = olympic.get_S_param(1, 1);


%% Calculate maximum gain in given direction

theta = 30*pi/180;
phi = 0*pi/180;
%olympic = olympic.set_source_max_rE_phi( theta, phi );
olympic = olympic.set_source_max_rE_theta( theta, phi );
antenna_group = 1;

% theta = 90*pi/180;
% phi = (90-45)*pi/180;
% olympic = olympic.set_source_max_rE_phi( theta, phi );
% %olympic = olympic.set_source_max_rE_theta( theta, phi );
% antenna_group = 2;

% theta = 90*pi/180;
% phi = 180*pi/180;
% olympic = olympic.set_source_max_rE_phi( theta, phi );
% %olympic = olympic.set_source_max_rE_theta( theta, phi );
% antenna_group = 3;

beam_pattern = olympic.array_beam_pattern(antenna_group);


% Plot Gain Pattern

aa = length(beam_pattern);
beam_pattern(aa).frequency;

figure
olympic.plot_surf_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_total) );
title(['realised gain theta at frequency ', num2str(frequency(aa)*1e-9), ' GHz']);

figure
%olympic.plot_spherical_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_total), -10 );
olympic.plot_spherical_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_total));


% % plot gain in particular direction vs frequency
% 
% [ realised_gain, frequency ] = HFSS_Tools.get_gain_total(beam_pattern, theta, phi);
% figure
% plot(frequency, olympic.dBp(realised_gain));
% xlabel('Frequency (Hz)');
% ylabel('Realised Gain (dBi)');
% title('get gain total');
% 
% 
% % plot gain in particular direction vs frequency
% 
% [ realised_gain, frequency ] = HFSS_Tools.get_gain_theta(beam_pattern, theta, phi);
% figure
% plot(frequency, olympic.dBp(realised_gain));
% xlabel('Frequency (Hz)');
% ylabel('Realised Gain (dBi)');
% title('get gain theta');
% 
% 
% % plot gain in particular direction vs frequency
% 
% [ realised_gain, frequency ] = HFSS_Tools.get_gain_phi(beam_pattern, theta, phi);
% figure
% plot(frequency, olympic.dBp(realised_gain));
% xlabel('Frequency (Hz)');
% ylabel('Realised Gain (dBi)');
% title('get gain phi');
% 
% 
% % plot maximum gain
% 
% [ realized_gain_phi, theta, phi, frequency ] = HFSS_Tools.get_max_gain_total(beam_pattern);
% figure
% plot(frequency, olympic.dBp(realized_gain_phi));
% xlabel('Frequency (Hz)');
% ylabel('Realised Gain (dBi)');
% title('get max gain total');
% 
% [ realized_gain_theta, theta, phi, frequency ] = HFSS_Tools.get_max_gain_theta(beam_pattern);
% figure
% plot(frequency, olympic.dBp(realized_gain_theta));
% xlabel('Frequency (Hz)');
% ylabel('E-Theta Gain (dBi)');
% title('get max gain theta');
% 
% [ realized_gain_phi, theta, phi, frequency ] = HFSS_Tools.get_max_gain_phi(beam_pattern);
% figure
% plot(frequency, olympic.dBp(realized_gain_phi));
% xlabel('Frequency (Hz)');
% ylabel('E-Phi Gain (dBi)');
% title('get max gain phi');
% 
% 
% Calculate beam width in theta and phi axis

[ beam_width_theta, beam_width_phi, frequency ] = HFSS_Tools.get_beam_width(beam_pattern);
figure
hold on
plot(frequency, beam_width_theta/HFSS_Tools.degrees,'r');
plot(frequency, beam_width_phi/HFSS_Tools.degrees,'b');
xlabel('Frequency (Hz)');
ylabel('Beam width (degrees)');
title('get beam width');

% Calculate the maximum sidelobe level

[ side_lobe_gain, theta_side, phi_side, frequency ] = HFSS_Tools.get_side_lobe(beam_pattern);
[ realised_gain, theta, phi, frequency ] = HFSS_Tools.get_max_gain_total(beam_pattern);
figure
hold on
plot(frequency, HFSS_Tools.dBp(side_lobe_gain),'r');
plot(frequency, HFSS_Tools.dBp(realised_gain),'b');
xlabel('Frequency (Hz)');
ylabel('Side lobe gain (dBi)');
title('get side lobe');


% plot the active S-parameters

[reflection_coefficient, frequency ] = olympic.array_active_Snn( antenna_group );
figure
hold on
plot(frequency, olympic.dBv(reflection_coefficient));
xlabel('Frequency (Hz)');
ylabel('Active Snn (dB)');
title('array active Snn');
grid on;

break;

% close all
%%

%%
figure

olympic.plot_surf_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_theta) );
title(['realised gain theta at frequency ', num2str(frequency(aa)*1e-9), ' GHz']);

figure
olympic.plot_spherical_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_theta), -30 );


[ beam_pattern2 ] = HFSS_Tools.rotate_axis(beam_pattern, 90*HFSS_Tools.degrees, 0, 0);

figure
olympic.plot_surf_pattern( beam_pattern2(aa).phi, beam_pattern2(aa).theta, HFSS_Tools.dBp(beam_pattern2(aa).realized_gain_phi) );
title(['realised gain theta at frequency ', num2str(frequency(aa)*1e-9), ' GHz']);

figure
olympic.plot_spherical_pattern( beam_pattern2(aa).phi, beam_pattern2(aa).theta, HFSS_Tools.dBp(beam_pattern2(aa).realized_gain_phi), -30 );
