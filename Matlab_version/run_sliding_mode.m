% AE 8803 VAM - Optimization-based Learning Control
% Project: RPO and docking using Adaptive Optimal Controller

global Jx Jy Jz Jtx Jty Jtz LcHIST

%% Input 

% time
t_max = 5e4;
timestep=1;
t_sim = 1:timestep:t_max;
nbrsteps = size(t_sim,2);

% Input motion of the asteroid
spin_axis = [1;0.5;1]; % intertiel
spin_axis = spin_axis/norm(spin_axis);
omega_rpm = 0.2;
omega_rads = omega_rpm/60*(2*pi);


% Initial attitude
e = cross(spin_axis,[1;0;0]); %rotation axis to define the quaternion
                % we assume that the capture system is along the x-axis on the sc body frame   
nu = 0;
q_c_init = [0; 0;0; 1]; % should orient the Capture system to the asteroid spin axis 
% (so its initial attitude depend on the spin axis of the asteroid but also on where it is on the orbit)
q_c_init = q_c_init/norm(q_c_init);

q_t_init = [0.0; 0; 0; 1];
q_t_init= q_t_init/norm(q_t_init);


% Quaternion inverse
q_t_inv_init = [-q_t_init(1:3);
                q_t_init(4)];

% Angular velocities
omega_c = [0; 0.0; 0.0];
omega_t = spin_axis*omega_rads; % Input: angular velocity

% Inertia matrix
 % assume SC has a shape of cylinder with the symmetry along x
R = 1.5;
H = 6;
m = 1e4; %(kg)
Jx = 1/2*m*R^2;  % (kg.m2)
Jy = 1/4*m*R^2 + 1/12*m*H^2;
Jz = 1/4*m*R^2 + 1/12*m*H^2;

 %assume the asteroid is a perfect sphere
mt = 1e6; % masse of the target (kg)
Rt = 5; % masse of the target (kg)
Jtx = 2/5*mt*R^2;
Jty = 2/5*mt*R^2;
Jtz = 2/5*mt*R^2;

% Propulsion charcateristics
Isp = 200;
thrust_nom = 5; % nominal thrust (N)
torque_nom_x = 1; % norminal torque for x-axis (Nm)
torque_nom_yz = 3; % norminal torque for y-axis, z-axis (Nm)

%% Program
% Initial quaternion error
delta_q_init = multiply_quat(q_c_init, q_t_inv_init);


y0 = [q_c_init; omega_c; q_t_init; omega_t];
skew_omega = cross_skew_matrix(omega_c);

%y_sim = ode4n_sliding_mode(@eulerseqns_attitude, t_sim, y0');
y_sim = ode4n_sliding_mode(@eulerseqns_attitude, t_sim, y0');

LxcHIST = LcHIST(1,:);
LycHIST = LcHIST(2,:);
LzcHIST = LcHIST(3,:);

% PWPF
K_p_x = 1/torque_nom_x; % proportional tuning gain for pwpf 
% normalize the command with respect to the nominal thrust of RCS
K_m_x = 20; % tuning gain for pwpf
T_m_x = 100; %tuning gain for pwpf (s)
U_on_x = 0.45; %tuning parameter for pwpf (schmitt trigger)
U_off_x = U_on_x*0.8; %tuning parameter for pwpf (schmitt trigger)

K_p_yz = 1/torque_nom_yz; % proportional tuning gain for pwpf 
% normalize the command with respect to the nominal thrust of RCS
K_m_yz = 20; % tuning gain for pwpf
T_m_yz = 100; %tuning gain for pwpf (s)
U_on_yz = 0.45;%0.9; %tuning parameter for pwpf (schmitt trigger)
U_off_yz = U_on_yz*0.8; %tuning parameter for pwpf (schmitt trigger)

% change the timestep for pwpf
% timestep_pwpf = timestep*0.1;
% t_sim_pwpf = 1:timestep_pwpf:t_max;
% nbrsteps_pwpf = size(t_sim_pwpf,2);

DC_all = zeros(3,1);
f_o_all = zeros(3,1);

[LxoHIST,DC_all(1,1),f_o(1,1)] = PWPF_Run(LxcHIST,K_p_x,K_m_x,T_m_x,U_on_x,U_off_x,t_sim);
[LyoHIST,DC(1,1),f_o(1,1)] = PWPF_Run(LycHIST,K_p_yz,K_m_yz,T_m_yz,U_on_yz,U_off_yz,t_sim);
[LzoHIST,DC(1,1),f_o(1,1)] = PWPF_Run(LzcHIST,K_p_yz,K_m_yz,T_m_yz,U_on_yz,U_off_yz,t_sim);

% verify result
%LHIST = M*UHIST*thrust_nom;
LoHIST = [LxoHIST*torque_nom_x;LyoHIST*torque_nom_yz;LzoHIST*torque_nom_yz];
y_sim_verif = ode4n_WithGivenTorqueControl(@eulerseqns_attitude, t_sim, y0', LoHIST);

% Relate the thrust and the torque, depending on the RCS configuration

% For the RCS, we assume the configuration 14 of the paper "A Study of
% Spacecraft Reaction Thruster Configurations  for Attitude Control System"
% We also assume that the center of gravity of the sc is on the plane formed by the ring of the RCS.
% Thus, L = CoM(x) in the matrix.

% In the following method, the torque are related to the thrust of each
% thruster, using the pseudo-inverse of the RCS configuration matrix
% the method used in the code can be applied for configuration such that for one thruster there is always another one that
% can produce thrust in the opposite direction from the same point.
% However, with more 'tricks', I think the method can be more generalized.
% The issue without this assumption is that for certain torque a thruster
% would have to fire in a negative direction or if not, the thursters at
% the opposite direction would fire alone, making the spacecraft translate.

%  RCS Configuration 14 in ' a study of sc RCS'
R = R; % R (m), is the thurters' arm radisu and is equal the radius of the spacecraft
L = 0; % L, the thrusters' arm lenght (m)
M = [ R, -R, 0, 0, R, -R, 0, 0, -R, R, 0, 0, -R, R, 0, 0;
      0,  0, -R, R, L, -L, 0, 0,  0, 0, R, -R, L, -L, 0, 0;
     -L, L,  0, 0, 0, 0, -R, R, -L, L, 0, 0, 0, 0, R, -R];
M_modified = [ R, 0, R, 0, -R, 0, -R, 0;
               0, -R, L, 0, 0, R, L, 0;
               -L, 0, 0, -R, -L, 0, 0, R]; % I deleted the columns of even number so that the thrust command can be negative

M_pinv = pinv(M_modified); %pseudo-inverse of M_modified

% the command and the output of the thrusts required to produce the torque along each axis
% ---------- Case without PWPF
TcHIST = M_pinv * LcHIST; 

% ---------- Case with PWPF
ToHIST = M_pinv * LoHIST; 
% We will not modify the thrust command and output TcHIST ToHIST such that
% the thrust commands are positive as we are only interested by the
% propellant mass

% prepare output
chaserAtt = zeros(length(y_sim),4);
targetAtt = zeros(length(y_sim),4);
chaserAngVel = zeros(length(y_sim),3);
targetAngVel = zeros(length(y_sim),3);
for i=1:nbrsteps
    chaserAtt(i,1:4) = y_sim(i,1:4);
    chaserAngVel(i,1:3) = y_sim(i,5:7);
    targetAtt(i,1) = y_sim(i,11);
    targetAtt(i,2:4) = y_sim(i,8:10);
    targetAngVel(i,1:3) = y_sim(i,12:14);
end

% Compute necessary time to match the spin rate
threshold = 1e-3;
T_f = -1; % time necessary to match the spin rate
counter1 = 0;
counter2 = 0;
for i=1:nbrsteps
    
    %fprintf('Chaser %.2f\n', chaserAngVel(i,:));
    %fprintf('Target %.2f\n', targetAngVel(i,:));
    %fprintf('difference %.2f\n', chaserAngVel(i,:)-targetAngVel(i,:));
    for j=1:3
        if abs(chaserAngVel(i,j)-targetAngVel(i,j)) < threshold
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
        T_f = i;
        break;
    end
end


%% Output --------------

% Calculate propellant consumption
g0 = 9.81;
TcHIST_abs = abs(TcHIST);
ToHIST_abs = abs(ToHIST);
mass_prop_fromCOMMAND = sum(TcHIST_abs(:))*timestep /(Isp*g0); 
mass_prop_fromOUTPUT = sum(ToHIST_abs(:))*timestep /(Isp*g0); %kg
disp('propellant mass consumption (kg)')
fprintf('Using the thrusters command (kg) %.2f\n', mass_prop_fromCOMMAND);
fprintf('Using the thrusters output (kg) %.2f\n', mass_prop_fromOUTPUT);
if T_f == -1
   disp('the chaser could not match the spin rate of the target')
else
    fprintf('Time necessary to match the spin rate (s) %.2f\n', T_f);
end

% Difference between control and output
disp('Tc - To at the end ')
diff = norm(LcHIST(:,end-1000) - LoHIST(:,end-1000))

% Maximum thrust
disp('maximum thrust')
thrust_max = max(TcHIST_abs)

save attitude.mat chaserAtt targetAtt

%% Plot ------------
% check if 
% figure()
% hold on
% plot(t_sim, chaserAngVel(:,1)-targetAngVel(:,1))
% plot(t_sim, chaserAngVel(:,2)-targetAngVel(:,2))
% plot(t_sim, chaserAngVel(:,3)-targetAngVel(:,3))
% ylabel('\omega_{c,i}-\omega_{t,i}')
% xlabel('Time [s]')
% title('The difference between the chaser (command) and target Attitude rate over visible time')
% legend('\omega_{c,1}-\omega_{c,1}','\omega_{c,2}-\omega_{c,2}','\omega_{c,3}-\omega_{c,3}','Location','best')
% ylim([-0.1,0.1])
% 
% figure()
% hold on
% plot(t_sim, y_sim(:,5))
% plot(t_sim, y_sim(:,6))
% plot(t_sim, y_sim(:,7))
% ylabel('\omega_{c,i}')
% xlabel('Time [s]')
% title('Chaser Attitude rate (command) over visible time')
% legend('\omega_{c,1}','\omega_{c,2}','\omega_{c,3}','Location','best')
% ylim([-0.1,0.1])
% 
% figure()
% hold on
% plot(t_sim, y_sim_verif(:,5))
% plot(t_sim, y_sim_verif(:,6))
% plot(t_sim, y_sim_verif(:,7))
% ylabel('\omega_{c,i}')
% xlabel('Time [s]')
% title('Chaser Attitude rate (output) over visible time')
% legend('\omega_{o,1}','\omega_{o,2}','\omega_{o,3}','Location','best')
% ylim([-0.1,0.1])
% 
% % figure();
% % hold on;
% % plot(t_sim, y_sim(:,5));
% % plot(t_sim, y_sim(:,6));
% % plot(t_sim(:,7));
% % plot(t_sim_pwpf, y_sim_verif(:,5));
% % plot(t_sim_pwpf, y_sim_verif(:,6));
% % plot(t_sim_pwpf, y_sim_verif(:,7));
% % ylabel('\omega_{c,i}');
% % xlabel('Time [s]');
% % title('Comparison between command and output of the chaser Attitude rate over visible time');
% % legend('\omega_{c,1}','\omega_{c,2}','\omega_{c,3}','\omega_{o,1}','\omega_{o,2}','\omega_{o,3}','Location','best');
% % ylim([-0.1,0.1]);
% 
% figure()
% hold on
% plot(t_sim, LcHIST(1,:))
% plot(t_sim, LcHIST(2,:))
% plot(t_sim, LcHIST(3,:))
% %plot(t_sim, sqrt(LcHIST(1,:).^2+LcHIST(2,:).^2+LcHIST(3,:).^2))
% ylabel('torque command (N.m)')
% xlabel('Time [s]')
% title('Torque command over visible time')
% legend('L_{c,1}','L_{c,2}','L_{c,3}','L', 'Location','best')
% %ylim([-1,1])
% 
% figure()
% hold on
% plot(t_sim, LoHIST(1,:))
% plot(t_sim, LoHIST(2,:))
% plot(t_sim, LoHIST(3,:))
% %plot(t_sim_pwpf, sqrt(LHIST(1,:).^2+LHIST(2,:).^2+LHIST(3,:).^2))
% ylabel('torque (N.m)')
% xlabel('Time [s]')
% title('Torque output over visible time')
% ylim([-torque_nom_x*1.1 torque_nom_x*1.1])
% legend('L_{c,1}','L_{c,2}','L_{c,3}','L', 'Location','best')

figure();

% First plot (Top Left)
subplot(2, 2, 1);
hold on;
plot(t_sim, y_sim(:,5));
plot(t_sim, y_sim(:,6));
plot(t_sim, y_sim(:,7));
ylabel('\omega_{c,i}');
xlabel('Time [s]');
title('Chaser Attitude rate (command) over visible time');
legend('\omega_{c,1}','\omega_{c,2}','\omega_{c,3}','Location','best');
ylim([-0.1, 0.1]);

% Second plot (Top Right)
subplot(2, 2, 2);
hold on;
plot(t_sim, y_sim_verif(:,5));
plot(t_sim, y_sim_verif(:,6));
plot(t_sim, y_sim_verif(:,7));
ylabel('\omega_{c,i}');
xlabel('Time [s]');
title('Chaser Attitude rate (output) over visible time');
legend('\omega_{o,1}','\omega_{o,2}','\omega_{o,3}','Location','best');
ylim([-0.1, 0.1]);

% Third plot (Bottom Left)
subplot(2, 2, 3);
hold on;
plot(t_sim, LcHIST(1,:));
plot(t_sim, LcHIST(2,:));
plot(t_sim, LcHIST(3,:));
%plot(t_sim, sqrt(LcHIST(1,:).^2+LcHIST(2,:).^2+LcHIST(3,:).^2));
ylabel('torque command (N.m)');
xlabel('Time [s]');
title('Torque command over visible time');
legend('L_{c,1}','L_{c,2}','L_{c,3}', 'Location','best');
ylim([-1, 1]);

% Fourth plot (Bottom Right)
subplot(2, 2, 4);
hold on;
plot(t_sim, LoHIST(1,:));
plot(t_sim, LoHIST(2,:));
plot(t_sim, LoHIST(3,:));
%plot(t_sim_pwpf, sqrt(LHIST(1,:).^2+LHIST(2,:).^2+LHIST(3,:).^2));
ylabel('torque (N.m)');
xlabel('Time [s]');
title('Torque output over visible time');
ylim([-torque_nom_x * 1.1, torque_nom_x * 1.1]);
legend('L_{o,1}','L_{o,2}','L_{o,3}', 'Location','best');

%figure()
%hold on
%plot(t_sim, u(:))
%ylabel('PWPF thruster On/Off (N.m)')
%xlabel('Time [s]')
%title('PWPF thruster On/Off over visible time')
%legend('u', 'Location','best')

% figure()
% hold on
% plot(t_sim, UcHIST(1,:))
% plot(t_sim, UcHIST(3,:))
% plot(t_sim, UcHIST(5,:))
% ylabel('thruter 1,3,5 commands (N)')
% xlabel('Time [s]')
% title('thruter 1,3,5 commands over time')
% legend('u_{1}','u_{2}','u_{3}', 'Location','best')
% 
% figure()
% hold on
% plot(t_sim_pwpf, UHIST(1,:))
% plot(t_sim_pwpf, UHIST(2,:))
% plot(t_sim_pwpf, UHIST(3,:))
% ylabel('thruter 1,2,3 control (N)')
% xlabel('Time [s]')
% title('thruter 1,2,3 output over time')
% legend('u_{1}','u_{2}','u_{3}', 'Location','best')
% 
% 
% 
for i = 2:3
    figure(); % Create a new figure window for each thruster
    hold on;
    plot(t_sim, LcHIST(i,:), 'b', 'LineWidth', 1);
    plot(t_sim, LoHIST(i,:), 'r', 'LineWidth', 1);
    ylabel(['Torque ' num2str(i) ' Torque (Nm)']);
    xlabel('Time [s]');
    title(['Torque ' num2str(i) ' Torque Command vs. On-Off Output']);
    legend('Continuous Command', 'On-Off Output', 'Location', 'best');

    grid on;
    hold off;
end






