clear all
close all

testproject = HFSS_Tools('Project8_forLen', 'forAnna', 12);
testproject = testproject.set_solution_setup('HFSS Setup 1');
testproject = testproject.set_solution_sweep_freq('Sweep 1');
testproject = testproject.set_solution_sweep_fields('Sweep 2');
testproject = testproject.set_source_type('Designer');
testproject = testproject.set_angle_resolution(1,1);

testproject = testproject.set_source_name( 1, 'port1');
testproject = testproject.set_source_name( 2, 'port2');
testproject = testproject.set_source_name( 3, 'port3');
testproject = testproject.set_source_name( 4, 'port4');
testproject = testproject.set_source_name( 5, 'port5');
testproject = testproject.set_source_name( 6, 'port6');
testproject = testproject.set_source_name( 7, 'port7');
testproject = testproject.set_source_name( 8, 'port8');
testproject = testproject.set_source_name( 9, 'port9');
testproject = testproject.set_source_name( 10, 'port10');
testproject = testproject.set_source_name( 11, 'port11');
testproject = testproject.set_source_name( 12, 'port12');

% % testproject = testproject.simulate_freq;
% % testproject = testproject.simulate_fields;

%% S-Parameters

testproject = testproject.DESIGNER_export_S_parameters;
testproject = testproject.import_S_parameters;


%% Fields

testproject = testproject.HFSS_export_fields_data;
testproject = testproject.import_fields_data;

%%

save('testprojectDesigner', 'testproject');

%%

% load('testprojectDesigner');

%%

testproject = testproject.set_source_groups([1,2,3,4,5,6,7,8], 1);
testproject = testproject.set_source_groups([9, 10, 11, 12], 2);

testproject.print_S_param(1);

% Plot S parameters

[frequency, S] = testproject.get_S_param(1, 1);


%% Calculate maximum gain in given direction

theta = 90*pi/180;
phi = 0*pi/180;
testproject = testproject.set_source_max_rE_phi( theta, phi );

antenna_group = 1;
beam_pattern = testproject.array_beam_pattern(antenna_group);


% Plot Gain Pattern

aa = length(beam_pattern);
beam_pattern(aa).frequency;

figure
testproject.plot_surf_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_total) );
title(['realised gain theta at frequency ', num2str(frequency(aa)*1e-9), ' GHz']);

figure
testproject.plot_spherical_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_total), -10 );


% plot gain in particular direction vs frequency

[ realised_gain, frequency ] = HFSS_Tools.get_gain_total(beam_pattern, theta, phi);
figure
plot(frequency, testproject.dBp(realised_gain));
xlabel('Frequency (Hz)');
ylabel('Realised Gain (dBi)');
title('get gain total');


% plot gain in particular direction vs frequency

[ realised_gain, frequency ] = HFSS_Tools.get_gain_theta(beam_pattern, theta, phi);
figure
plot(frequency, testproject.dBp(realised_gain));
xlabel('Frequency (Hz)');
ylabel('Realised Gain (dBi)');
title('get gain theta');


% plot gain in particular direction vs frequency

[ realised_gain, frequency ] = HFSS_Tools.get_gain_phi(beam_pattern, theta, phi);
figure
plot(frequency, testproject.dBp(realised_gain));
xlabel('Frequency (Hz)');
ylabel('Realised Gain (dBi)');
title('get gain phi');


% plot maximum gain

[ realized_gain_phi, theta, phi, frequency ] = HFSS_Tools.get_max_gain_total(beam_pattern);
figure
plot(frequency, testproject.dBp(realized_gain_phi));
xlabel('Frequency (Hz)');
ylabel('Realised Gain (dBi)');
title('get max gain total');

[ realized_gain_theta, theta, phi, frequency ] = HFSS_Tools.get_max_gain_theta(beam_pattern);
figure
plot(frequency, testproject.dBp(realized_gain_theta));
xlabel('Frequency (Hz)');
ylabel('E-Theta Gain (dBi)');
title('get max gain theta');

[ realized_gain_phi, theta, phi, frequency ] = HFSS_Tools.get_max_gain_phi(beam_pattern);
figure
plot(frequency, testproject.dBp(realized_gain_phi));
xlabel('Frequency (Hz)');
ylabel('E-Phi Gain (dBi)');
title('get max gain phi');


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
figure
hold on
plot(frequency, theta,'r');
plot(frequency, phi,'b');
plot(frequency, theta_side,'b');
plot(frequency, phi_side,'g');
xlabel('Frequency (Hz)');
ylabel('Side lobe gain (dBi)');
title('get side lobe');


% plot the active S-parameters

[reflection_coefficient, frequency ] = testproject.array_active_Snn( antenna_group );
figure
hold on
plot(frequency, testproject.dBv(reflection_coefficient));
xlabel('Frequency (Hz)');
ylabel('Active Snn (dB)');
title('array active Snn');

% close all
%%

%%
figure

testproject.plot_surf_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_theta) );
title(['realised gain theta at frequency ', num2str(frequency(aa)*1e-9), ' GHz']);

figure
testproject.plot_spherical_pattern( beam_pattern(aa).phi, beam_pattern(aa).theta, HFSS_Tools.dBp(beam_pattern(aa).realized_gain_theta), -30 );


[ beam_pattern2 ] = HFSS_Tools.rotate_axis(beam_pattern, 90*HFSS_Tools.degrees, 0, 0);

figure
testproject.plot_surf_pattern( beam_pattern2(aa).phi, beam_pattern2(aa).theta, HFSS_Tools.dBp(beam_pattern2(aa).realized_gain_phi) );
title(['realised gain theta at frequency ', num2str(frequency(aa)*1e-9), ' GHz']);

figure
testproject.plot_spherical_pattern( beam_pattern2(aa).phi, beam_pattern2(aa).theta, HFSS_Tools.dBp(beam_pattern2(aa).realized_gain_phi), -30 );
