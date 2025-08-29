% for VGC ACM, RPO phase, Step 3 (match the asteroid spin rate)
%  by Akira Jinkoji
% Based on the code of Alexandre Masset




%% Map the propellant mass and Tof with respect to the spin axis 
global resultsHIST
% Output
    % max_prop_mass: maximum propellant mass consumed with respect to the spin
    % axis of the asteroid, given a constant spin rate
    % case_maxMass: [spin axis, spin rate, mass_prop, ToF], the input and the
    % outputs for the case with the max propellant consumption with respect to the spin
    % axis of the asteroid, given a constant spin rate  
    % max_tof, case_maxTof: same principle but for the maximum Time of
    % flight (ToF)

% Input
mass_t = 1.45e6; % mass of the asteroid / target (kg)
SpinRate_t = 0.2; % spin rate of the asteroid / target (rpm)

mass_sc = 10244; %mass of the spacecraft 

x_values = -1:0.4:1; %maybe to be modified because symmetric
y_values = -1:0.4:1;
z_values = -1:0.4:1;
nbrsteps = numel(x_values)*numel(y_values)*numel(z_values);

% initialize storage of results
resultsHIST = zeros(6,nbrsteps);
counter = 0;

%SpinAxis_t = [1;0;0];
%[prop_mass, Tof] = match_spin_rate(SpinRate_t, SpinAxis_t, mass_t, mass_sc);
for x = x_values
    for y = y_values
        for z = z_values
            if x ~= 0 || y ~= 0 || z ~= 0 % Exclude the origin
                counter = counter + 1;
                fprintf('Simulation %.2f\n', counter);
                SpinAxis_t = [x;y;z]/norm([x;y;z]);
                [prop_mass, Tof, thurst_max] = match_spin_rate_step3(SpinRate_t, SpinAxis_t, mass_t, mass_sc);
                resultsHIST(1:3, counter) = SpinAxis_t;
                resultsHIST(4, counter) = prop_mass;
                resultsHIST(5, counter) = Tof;
                resultsHIST(6,counter) = thurst_max;
            end
        end
    end
end

% Remove unused preallocated space
resultsHIST = resultsHIST(:, 1:counter);

% Finding the maximum propellant mass and corresponding index
[max_prop_mass, index_maxProp] = max(resultsHIST(4,:));
case_maxMass = [resultsHIST(:,index_maxProp)]

% Finding the maximum time of flight and corresponding index
[max_ToF, index_maxToF] = max(resultsHIST(5,:));
case_maxToF = [resultsHIST(:,index_maxToF)]

% Convert to spherical coordinates
theta = atan2(resultsHIST(2, :), resultsHIST(1, :));
phi = zeros(1,size(resultsHIST,2));
for j = 1:size(resultsHIST,2)
    phi(j) = acos(resultsHIST(3) / sqrt(resultsHIST(1,j)^2 + resultsHIST(2,j)^2 + resultsHIST(3,j)^2));
end

% Create a full-sphere grid for theta and phi
theta_grid = linspace(-pi, pi, numel(theta)); % Full azimuth range
phi_grid = linspace(0, pi, numel(phi)); % Full elevation range

[Theta, Phi] = meshgrid(theta_grid, phi_grid);
% Interpolate propellant mass for the grid
F_mass = scatteredInterpolant(theta', phi', resultsHIST(4, :)','linear', 'none');
Z_mass = F_mass(Theta, Phi);

F_tof = scatteredInterpolant(theta', phi', resultsHIST(5, :)','linear', 'none');
Z_tof = F_tof(Theta, Phi);

% Convert spherical to Cartesian coordinates for plotting
[X, Y, Z_sphere] = sph2cart(Theta, Phi, 1);

% Plot the surface with interpolated values
figure;
surf(X, Y, Z_sphere, Z_mass, 'EdgeColor', 'none');
colorbar;
c = colorbar;
c.Label.String = 'Propellant Mass Consumed (kg)';
title('Propellant Mass Consumed vs Spin Axis Direction');
xlabel('X');
ylabel('Y');
zlabel('Z');
colormap jet;
shading interp;
axis equal;
% Add a text box with the input parameters
%annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', sprintf('mass of the asteroid = %.2e kg ; target spin rate = %.2f rpm', mass_t, SpinRate_t), 'FitBoxToText', 'on', 'BackgroundColor', 'white')

% Plot the surface with interpolated values
figure;
surf(X, Y, Z_sphere, Z_tof, 'EdgeColor', 'none');
colorbar;
c = colorbar;
c.Label.String = 'Time of flight (sec)';
title('Time of flight vs Spin Axis Direction');
xlabel('X');
ylabel('Y');
zlabel('Z');
colormap jet;
shading interp;
axis equal;
%annotation('textbox', [0.02, 0.6, 0.3, 0.25], 'String', sprintf('mass of the asteroid = %.2e kg ; target spin rate = %.2f rpm', mass_t, SpinRate_t), 'FitBoxToText', 'on', 'BackgroundColor', 'white')

% Plot the 3D scatter plot
figure;
scatter3(resultsHIST(1, :), resultsHIST(2, :), resultsHIST(3, :), 50, resultsHIST(4, :), 'filled');
colorbar;
c = colorbar;
c.Label.String = 'Propellant Mass Consumed (kg)';
title('Propellant mass consumed VS Spin axis direction');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis 

%% Propellant mass consumption + time of flight with respect to the SC mass
clear all 
global results2HIST

%Input
nbrsteps = 20;
mass_sc = linspace(1e3, 3e4,nbrsteps);
mass_t = 1e6;
SpinRate_t = 0.2; %rpm
SpinAxis_t = [0.5;0.5;0.7];


% initialize storage of results
results2HIST = zeros(4,nbrsteps);
counter = 0;

for mass = mass_sc
    counter = counter + 1
    [prop_mass, Tof,thrustmax] = match_spin_rate_step3(SpinRate_t, SpinAxis_t, mass_t, mass);
    results2HIST(1,counter) = mass;
    results2HIST(2, counter) = prop_mass;
    results2HIST(3, counter) = Tof;
    results2HIST(4, counter) = thrustmax;
end

figure()
hold on
plot(mass_sc, results2HIST(2,:))
ylabel('Consumed propellant mass (kg)')
xlabel('Spacecraft mass (kg)')
title('Propellant mass consumption VS Spacecraft mass')
%legend(', 'Location','best')

figure()
hold on
plot(mass_sc, results2HIST(4,:))
ylabel('maximum thrust command (N)')
xlabel('mass of the spacecraft (kg)')
title('Mass of the asteroid VS maximum thrust command')
ylim([0,0.01])
%legend('', 'Location','best')

figure()
hold on
plot(mass_sc, results2HIST(3,:))
ylabel('Time of flight (sec)')
xlabel('Spacecraft mass (kg)')
title('Time of flight VS Spacecraft mass')
%legend(', 'Location','best')



%% Propellant mass consumption + time of flight with respect to the spin rate

clear all 
global results3HIST

%Input
nbrsteps = 50;
mass_sc = 1e3;
mass_t = 1e6;
SpinRate_t = linspace(0,0.5,nbrsteps); %rpm
SpinAxis_t = [0.5;0.5;0.7];

% initialize storage of results
results3HIST = zeros(3,nbrsteps);
counter = 0;

for spin_rate = SpinRate_t
    counter = counter + 1
    [prop_mass, Tof] = match_spin_rate_step3(spin_rate, SpinAxis_t, mass_t, mass_sc)
    results3HIST(1,counter) = spin_rate;
    results3HIST(2, counter) = prop_mass;
    results3HIST(3, counter) = Tof;
end

figure()
hold on
plot(SpinRate_t, results3HIST(2,:))
ylabel('Consumed propellant mass (kg)')
xlabel('Spin rate of the asteroid (rpm)')
title('Propellant mass consumption VS Asteroid spin rate')
%legend(', 'Location','best')

figure()
hold on
plot(SpinRate_t, results3HIST(3,:))
ylabel('Time of flight (kg)')
xlabel('Spin rate of the asteroid (rpm)')
title('Time of flight VS Asteroid spin rate')
%legend(', 'Location','best')

%% Test 

%Input
nbrsteps = 50;
mass_sc = 1e4;
mass_t = 1e6;
SpinRate_t = 0.2; %rpm
SpinAxis_t = [0.1;-0.5;-0.9];

[prop_mass, Tof, thrust_max] = match_spin_rate_step3(SpinRate_t, SpinAxis_t, mass_t, mass_sc);