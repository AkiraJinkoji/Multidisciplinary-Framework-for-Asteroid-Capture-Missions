% for VGC ACM, RPO phase, Step 3 (match the asteroid spin rate)
%  by Akira Jinkoji
% Based on the code of Alexandre Masset

% Input
% target_SpinRate : (rpm)
% target_SpinAxis: vector 3*1 in inertial frame 
% target_Mass: (kg)

% Output
    % prop_mass: the total propellant consumed for the operation (kg) 
    % Tof: total mass required for the operation (sec)
function [mass_prop, ToF,thrust_max] = match_spin_rate_step3(SpinRate_t, SpinAxis_t, mass_t, mass_sc)

global Jx Jy Jz Jtx Jty Jtz LcHIST

%% Input 

% time
t_max_reg = 1e5; % for regulation
timestep=1;
t_sim_reg = 1:timestep:t_max_reg;
nbrsteps_reg = size(t_sim_reg,2);

t_max_track = 5e4;
timestep=1;
t_sim_track = 1:timestep:t_max_track;
nbrsteps_track = size(t_sim_track,2);

t_max_total = t_max_track + t_max_reg;
t_sim_total = 1:timestep:t_max_total;
nbrsteps_total = size(t_sim_total,2);

% Input motion of the asteroid
spin_axis = SpinAxis_t; % inertiel
spin_axis = spin_axis/norm(spin_axis);
omega_rpm = SpinRate_t;
omega_rads = omega_rpm/60*(2*pi);


% Initial attitude
q_c_init = [0; 0;0; 1]; % neutral position
%x_axis_t = [1;0;0];
%q_c_init = compute_quaternion_attitude_from_2vec(-x_axis_t, spin_axis); 
% (so its initial attitude depend on the spin axis of the asteroid but also on where it is on the orbit)
q_c_init = q_c_init/norm(q_c_init);

%q_t_init = [0.0; 0; 0; 1];
x_axis_t = [1;0;0];
q_t_init = compute_quaternion_attitude_from_2vec(-x_axis_t, spin_axis); 
% attitude intertial --> body frame such that  it convert the assumed spin axis in inertial frame to minus the x-axis of asteroid in the body frame
q_t_init= q_t_init/norm(q_t_init);

spin_axis_BodyF = compute_vector_from_quaternion_attitude(q_t_init,spin_axis);


% Quaternion inverse
q_t_inv_init = [-q_t_init(1:3);
                q_t_init(4)];

% Angular velocities
omega_c = [0; 0.0; 0.0];
%omega_t = spin_axis*omega_rads; % Input: 
omega_t = [1;0;0]*omega_rads; % angular velocity, as we assume the x-axis of the asteroid in its body frame is its spin axis in inertial frame

% Inertia matrix
 % assume SC has a shape of cylinder with the symmetry along x
R = 1.5;
H = 6;
m = mass_sc; %(kg)
Jx = 1/2*m*R^2;  % (kg.m2)
Jy = 1/4*m*R^2 + 1/12*m*H^2;
Jz = 1/4*m*R^2 + 1/12*m*H^2;

 %assume the asteroid is a perfect sphere
mt = mass_t; % masse of the target (kg)
Rt = 5; % masse of the target (kg)
Jtx = 2/5*mt*R^2;
Jty = 2/5*mt*R^2;
Jtz = 2/5*mt*R^2;

% Propulsion charcateristics
Isp = 200;
max_acceptable_thrust = 100; %(N) 

% for pwpf
thrust_nom = 5; % nominal thrust (N)
torque_nom_x = 1; % norminal torque for x-axis (Nm)
torque_nom_yz = 10; % norminal torque for y-axis, z-axis (Nm)

% Max iter
max_iter = 5;

% initialize parameters for sliding mode
% Initial parameters for the regulation control
k_d = 1e-3; %gain for the spin rate to reach the sliding mode
k_p_reg =k_d/1000;   % gain for the attitude to reach the sliding surface (the bigger it is, the faster the convergence but the bigger the chattering)

% Initial parameters for the tracking control
k_p_track = 1e-3;   % gain to reach the sliding surface (the bigger it is, the faster the convergence but the bigger the chattering)
k_s = 1e-6;%7; % gain for the saturation function (the smaller it is, the smaller the chattering, limiting the max torque)  


%% Program

% iteration to limit the max_thrust
for iteration = 1:max_iter

% ------ For regulartion: bring the x-axis of the body frame of the spacecraft points to the spin axis of the asteroid 
% Initial quaternion error
y0_reg = [q_c_init; omega_c; q_t_init];
y_sim_reg = ode4n_sliding_mode_regulation(@eulerseqns_attitude_regulation, t_sim_reg, y0_reg', k_d, k_p_reg);

LxcHIST_reg = LcHIST(1,:);
LycHIST_reg = LcHIST(2,:);
LzcHIST_reg = LcHIST(3,:);

% prepare output
chaserAtt_reg = zeros(length(y_sim_reg),4);
targetAtt_reg = zeros(length(y_sim_reg),4);
chaserAngVel_reg = zeros(length(y_sim_reg),3);
for i=1:nbrsteps_reg
    chaserAtt_reg(i,1:4) = y_sim_reg(i,1:4);
    chaserAngVel_reg(i,1:3) = y_sim_reg(i,5:7);
    targetAtt_reg(i,1:4) = y_sim_reg(i,8:11);
end

% Compute necessary time to achieve regulation
threshold_Att = 1e-3;
T_f_reg = -1; % time necessary to match the spin rate
counter1 = 0;
counter2 = 0;
for i=1:nbrsteps_reg
    for j = 1:4
        if abs((chaserAtt_reg(i,j)-targetAtt_reg(i,j))) < threshold_Att
            counter1 = counter1 + 1;
        end
    end
    if counter1 == 4
        counter2 = counter2 + 1;
        counter1 = 0;
    else
        counter1 = 0;
        counter2 = 0;
    end
    
    if counter2 == 60
        T_f_reg = i;
        break;
    end
end
 
% ---- Tracking: match the spin rate of the asteroid once the spacecraft is
% oriented correctly with respect to the asteroid spin axis
q_c_init_track = chaserAtt_reg(end,:);
omega_c_track = chaserAngVel_reg(end, :);
q_t_init_track = q_t_init; % i assumed that the asteroid was not spinning during the regulation for simplicity
 
y0_track = [q_c_init_track'; omega_c_track'; q_t_init_track; omega_t];
y_sim_track = ode4n_sliding_mode_tracking(@eulerseqns_attitude, t_sim_track, y0_track', k_p_track, k_s);

LxcHIST_track = LcHIST(1,1:nbrsteps_track);
LycHIST_track = LcHIST(2,1:nbrsteps_track);
LzcHIST_track = LcHIST(3,1:nbrsteps_track);

LxcHIST = cat(2,LxcHIST_reg, LxcHIST_track);
LycHIST = cat(2, LycHIST_track, LycHIST_reg);
LzcHIST = cat(2, LzcHIST_reg, LzcHIST_track);

LcHIST_total = cat(1, LxcHIST, LycHIST);
LcHIST_total = cat(1, LcHIST_total, LzcHIST);

% prepare output
chaserAtt_track = zeros(length(y_sim_track),4);
targetAtt_track = zeros(length(y_sim_track),4);
chaserAngVel_track  = zeros(length(y_sim_track),3);
targetAngVel_track  = zeros(length(y_sim_track),3);
for i=1:nbrsteps_track
    chaserAtt_track(i,1:4) = y_sim_track(i,1:4);
    chaserAngVel_track(i,1:3) = y_sim_track(i,5:7);
    targetAtt_track(i,1:4) = y_sim_track(i,8:11);
    targetAngVel_track(i,1:3) = y_sim_track(i,12:14);
end

% Compute necessary time to achieve regulation
threshold_SpinRate = 1e-4;
T_f_track = -1; % time necessary to match the spin rate
counter1 = 0;
counter2 = 0;
counter3 = 0;
for i=1:nbrsteps_track
    for j=1:3
        if abs(chaserAngVel_track(i,j)-targetAngVel_track(i,j)) < threshold_SpinRate 
            counter1 = counter1 + 1;
        end
    end
    if counter1 == 3
        counter2 = counter2 + 1;
        counter1 = 0;
    else
        counter1 = 0;
        counter2 = 0;
    end
    
    if counter2 == 60
        T_f_track = i;
        break;
    end
end

% Concatenate the state during the regulation and tracking
chaserAngVel_total = cat(1, chaserAngVel_reg, chaserAngVel_track);
chaserAtt_total = cat(1, chaserAtt_reg,chaserAtt_track);
%targetAtt_total = cat(1, targetAtt_track(1:nbrsteps_reg,:),targetAtt_track);



% PWPF
K_p_x = 1/torque_nom_x; % proportional tuning gain for pwpf 
% normalize the command with respect to the nominal thrust of RCS
K_m_x = 9; % tuning gain for pwpf
T_m_x = 900; %tuning gain for pwpf (s)
U_on_x = 0.8;%0.9; %tuning parameter for pwpf (schmitt trigger)
U_off_x = U_on_x*0.3; %tuning parameter for pwpf (schmitt trigger)

K_p_yz = 1/torque_nom_yz; % proportional tuning gain for pwpf 
% normalize the command with respect to the nominal thrust of RCS
K_m_yz = 23; % tuning gain for pwpf
T_m_yz = 850; %tuning gain for pwpf (s)
U_on_yz = 0.7;%0.9; %tuning parameter for pwpf (schmitt trigger)
U_off_yz = U_on_yz*0.3; %tuning parameter for pwpf (schmitt trigger)

% change the timestep for pwpf
% timestep_pwpf = timestep*0.1;
% t_sim_pwpf = 1:timestep_pwpf:t_max;
% nbrsteps_pwpf = size(t_sim_pwpf,2);

DC_all = zeros(3,1);
f_o_all = zeros(3,1);

% [LxoHIST,DC_all(1,1),f_o(1,1)] = PWPF_Run(LxcHIST,K_p_x,K_m_x,T_m_x,U_on_x,U_off_x,t_sim_total);
% [LyoHIST,DC(1,1),f_o(1,1)] = PWPF_Run(LycHIST,K_p_yz,K_m_yz,T_m_yz,U_on_yz,U_off_yz,t_sim_total);
% [LzoHIST,DC(1,1),f_o(1,1)] = PWPF_Run(LzcHIST,K_p_yz,K_m_yz,T_m_yz,U_on_yz,U_off_yz,t_sim_total);
% 
% % verify result
% %LHIST = M*UHIST*thrust_nom;
% LoHIST = [LxoHIST*torque_nom_x;LyoHIST*torque_nom_yz;LzoHIST*torque_nom_yz];
% y_sim_verif = ode4n_WithGivenTorqueControl(@eulerseqns_attitude, t_sim, y0', LoHIST);

% Relate the thrust and the torque, depending on the RCS configuration

% For the RCS, we assume the configuration 14 of the paper "A Study of
% Spacecraft Reaction Thruster Configurations  for Attitude Control System"
% We also assume that the center of gravity of the sc is on the plane formed by the ring of the RCS.
% Thus, L = CoM(x) in the matrix.
% 
% In the following method, the torque are related to the thrust of each
% thruster, using the pseudo-inverse of the RCS configuration matrix
% the method used in the code can be applied for configuration such that for one thruster there is always another one that
% can produce thrust in the opposite direction from the same point.
% However, with more 'tricks', I think the method can be more generalized.
% The issue without this assumption is that for certain torque a thruster
% would have to fire in a negative direction or if not, the thursters at
% the opposite direction would fire alone, making the spacecraft translate.

%  RCS Configuration 14 in ' a study of sc RCS'
R = R; % R, is the thurters' arm radisu and is equal the radius of the spacecraft
L = 0 ; % L, the thrusters' arm lenght 
M = [ R, -R, 0, 0, R, -R, 0, 0, -R, R, 0, 0, -R, R, 0, 0;
      0,  0, -R, R, L, -L, 0, 0,  0, 0, R, -R, L, -L, 0, 0;
     -L, L,  0, 0, 0, 0, -R, R, -L, L, 0, 0, 0, 0, R, -R];
M_modified = [ R, 0, R, 0, -R, 0, -R, 0;
               0, -R, L, 0, 0, R, L, 0;
               -L, 0, 0, -R, -L, 0, 0, R]; % I deleted the columns of even number so that the thrust command can be negative

M_pinv = pinv(M_modified); %pseudo-inverse of M_modified

% the command and the output of the thrusts required to produce the torque along each axis
% ---------- Case without PWPF
TcHIST = M_pinv * LcHIST_total; 

% ---------- Case with PWPF
%ToHIST = M_pinv * LoHIST; 
% We will not modify the thrust command and output TcHIST ToHIST such that
% the thrust commands are positive as we are only interested by the
% propellant mass


%% Output --------------
% Calculate propellant consumption
g0 = 9.81;
TcHIST_abs = abs(TcHIST);
%ToHIST_abs = abs(ToHIST);
mass_prop_fromCOMMAND = sum(TcHIST_abs(:))*timestep /(Isp*g0); 
%mass_prop_fromOUTPUT = sum(ToHIST_abs(:))*timestep /(Isp*g0); %kg

mass_prop = mass_prop_fromCOMMAND;
ToF = t_max_reg + T_f_track;

% Difference between control and output
%disp('Tc - To at the end ')
%diff = norm(LcHIST(:,end-1000) - LoHIST(:,end-1000))

% Maximum thrust
thrust_max_reg = 0;
thrust_max_i = 0;
for i = size(TcHIST_abs,1)
    thrust_max_i = max(TcHIST_abs(i,1:size(LxcHIST_reg,2)));
    if thrust_max_i > thrust_max_reg
        thrust_max_reg = thrust_max_i;
    end
end
thrust_max_track = 0;
thrust_max_i = 0;
for i = size(TcHIST_abs,1)
    thrust_max_i = max(TcHIST_abs(i,size(LxcHIST_reg,2)+1:end));
    if thrust_max_i > thrust_max_track
        thrust_max_track = thrust_max_i;
    end
end
thrust_max = max([thrust_max_track, thrust_max_reg]);

% Condition to exit for loop: check if max_thrust < max_acceptable_thrust
if thrust_max < max_acceptable_thrust && T_f_track > -1 && T_f_reg > -1
    fprintf('Iteration %.2f\n', iteration);
    disp('Spacecraft matched the spin rate and respected the max thrust constraint')
    fprintf('Propellant mass consumed %.2f\n', mass_prop);
    %fprintf('DeltaV %.2f\n', iteration);
    fprintf('Time of flight %.2f\n', ToF);
    fprintf('Maximum thrust %.2f\n', thrust_max);
    break
else
    fprintf('Iteration %.2f\n', iteration);
    if thrust_max_reg >= max_acceptable_thrust
        disp('condition for max thrust not respected')
        k_p_reg = k_p_reg/10;
        k_d = k_d/10;
    end
    if thrust_max_track >= max_acceptable_thrust
        disp('condition for max thrust not respected')
        k_p_track = k_p_track/10;
        
    end
    if T_f_reg == -1 
        disp('spacecraft could not achieve the regulation within the given time')
        % time
        t_max_reg = t_max_reg*5;
        timestep=1;
        t_sim_reg = 1:timestep:t_max_reg;
        nbrsteps_reg = size(t_sim_reg,2);
        t_max_total = t_max_track + t_max_reg;
        t_sim_total = 1:timestep:t_max_total;
        nbrsteps_total = size(t_sim_total,2);
    end
    if T_f_track == -1 
        disp('spacecraft could not achieve the tracking within the given time')
        % time
        t_max_track = t_max_track*5;
        timestep=1;
        t_sim_track = 1:timestep:t_max_track;
        nbrsteps_track = size(t_sim_track,2);
        t_max_total = t_max_track + t_max_reg;
        t_sim_total = 1:timestep:t_max_total;
        nbrsteps_total = size(t_sim_total,2);
    end
end

if iteration == max_iter
    error('the control law could not match the spin rate within the given time and the max thrust constraint')
end 

end

%% verify the rotation axis of target and sc

ast_spin_axis_in_body_frame = zeros(3,nbrsteps_total);
for i = 1:nbrsteps_total
    % qt = targetAtt_track(i,:);
    % qt = qt/norm(qt); % normalize
    qc = chaserAtt_total(i,:);
    qc = qc/norm(qc); % normalize
    % Rotate v2 using the quaternion
    %SpinAxis_4t_BodyF = compute_vector_from_quaternion_attitude(qt,spin_axis/norm(spin_axis));
    SpinAxis_4c_BodyF = compute_vector_from_quaternion_attitude(qc,spin_axis/norm(spin_axis));
    % Extract the vector part of the result
    %target_spin_axis_in_body_frame(1:3,i) = SpinAxis_4t_BodyF;
    ast_spin_axis_in_body_frame(1:3,i) = SpinAxis_4c_BodyF;
end

%% Plot

figure();
% 
% First plot (Top Left)
subplot(2, 2, 1);
hold on;
plot(t_sim_total, chaserAngVel_total(:,1));
plot(t_sim_total, chaserAngVel_total(:,2));
plot(t_sim_total, chaserAngVel_total(:,3));
ylabel('\omega_{c,i}');
xlabel('Time [s]');
title('Chaser Attitude rate (command) over visible time');
legend('\omega_{c,1}','\omega_{c,2}','\omega_{c,3}','Location','best');
ylim([-0.1, 0.1]);

% Second plot (Top Right)
subplot(2, 2, 2);
hold on;
plot(t_sim_total, ast_spin_axis_in_body_frame(1,:));
plot(t_sim_total, ast_spin_axis_in_body_frame(2,:));
plot(t_sim_total, ast_spin_axis_in_body_frame(3,:));
ylabel('Asteroid spin axis in the spacecraft body frame');
xlabel('Time [s]');
title('Check if the Capture System (x-axis) is pointing toward the asteroid spin axis');
legend('$\mathrm{SpinAxis}_{x}$', '$\mathrm{SpinAxis}_{y}$', '$\mathrm{SpinAxis}_{z}$', 'Location', 'best', 'Interpreter', 'latex')

% Third plot (Bottom Left)
subplot(2, 2, 3);
hold on;
plot(t_sim_total, LcHIST_total(1,:));
plot(t_sim_total, LcHIST_total(2,:));
plot(t_sim_total, LcHIST_total(3,:));
%plot(t_sim, sqrt(LcHIST(1,:).^2+LcHIST(2,:).^2+LcHIST(3,:).^2));
ylabel('torque command (N.m)');
xlabel('Time [s]');
title('Torque command over visible time');
legend('L_{c,1}','L_{c,2}','L_{c,3}', 'Location','best');
%ylim([-1, 1]);
% 
% Fourth plot
subplot(2, 2, 4);
hold on
plot(t_sim_total, chaserAtt_total(:,1))
plot(t_sim_total, chaserAtt_total(:,2))
plot(t_sim_total, chaserAtt_total(:,3))
plot(t_sim_total, chaserAtt_total(:,4))
% plot(t_sim_total, ones(nbrsteps_total, 1)*q_t_init(1))
% plot(t_sim_total, ones(nbrsteps_total, 1)*q_t_init(2))
ylabel('q_i')
xlabel('Time [s]')
title('Chaser Attitude quaternion over visible time')
legend('q_1','q_2','q_3','q_4','Location','best')
ylim([-1,1])

