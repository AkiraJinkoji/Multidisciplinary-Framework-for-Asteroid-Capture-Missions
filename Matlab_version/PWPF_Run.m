%pwpf matlab (from simulink design)
%model based off of Trond Dagfinn Krovel's Paper: "Optimal Selection of PWPF Parameters"

% Modified by Akira for its VGC
function [u,DC,f_o] = PWPF_Run(rcHIST,K_p,K_m,T_m,U_on,U_off,t_sim)

%initialize
f_zero = 0;
U_initial = 0;
b = 0;
f = [];
u = [];
uconcat = U_initial;
u = [u uconcat];

NbrPulse = 0;
counter = 1;
for t = t_sim
    e = K_p*rcHIST(1,counter)-u(counter); 
    f_time = f_zero + (K_m*e - f_zero)*(1 - exp(-(t-b)/T_m));

    
    if f_time > 0
         if f_time < U_off
             uconcat = 0;
             f_zero = U_off;
             f_time = U_off;
             b = t;
         elseif f_time > U_on
             uconcat = 1;
             f_zero = U_on;
             f_time = U_on;
             b = t;
         else
             uconcat = u(end); %keep the same control than the previous one
             %f_zero = f_time; %added by Akira
             %b =t;% added by Akira
         end
    else
        
         if f_time > -U_off
             uconcat = 0;
             f_zero = -U_off;
             f_time = -U_off;
             b = t;
         elseif f_time < -U_on
             uconcat = -1;
             f_zero = -U_on;
             f_time = -U_on;
             b = t;
         else
             uconcat = u(end); %keep the same control than the previous one
             %f_zero = f_time; %added by Akira
             %b =t;% added by Akira
         end
    end
     f =[f f_time];
     NbrPulse = NbrPulse + uconcat;
     u = [u uconcat];
     counter = counter + 1;
 end
u = u(1:end-1);
h = U_on - U_off;
Ton = -T_m*log(1- h/(U_on - K_m*(K_p*rcHIST(1,end)-u(end)))); % <---- probably wrong formula
Toff = -T_m*log(1-h/(K_m*K_p*rcHIST(1,end)-U_off));% <---- probably wrong formula
DC = Ton/(Ton + Toff);
f_o = 1/(Ton + Toff);

%% Plot 
% 
% figure()
% hold on
% plot(t_sim, f(1,:))
% %plot(t_sim, u(1,:))
% plot(t_sim, U_off *ones(size(t_sim)))
% plot(t_sim, U_on* ones(size(t_sim)))
% plot(t_sim, -U_off *ones(size(t_sim)))
% plot(t_sim, -U_on* ones(size(t_sim)))
% ylabel(' ')
% xlabel('Time [s]')
% title('pwpf over time')
% legend('f(t)','u','u_{off}','u_{on}' ,'Location','best')

end