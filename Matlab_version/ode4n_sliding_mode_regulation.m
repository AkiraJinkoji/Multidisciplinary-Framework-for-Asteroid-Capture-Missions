function Y = ode4n_sliding_mode_regulation(odefun,tspan,y0,k_d, k_p, varargin)
% ODE4N modified ODE4 fixed step size integrator
% 
% Most recently modified by Akira Jinkoji, for the step 6 
% of the RPO phase of VGC ACM
% Modified version of ode4n provided by Dr. Glenn Lightsey in by Alexandre
% Masset

%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates 
%   the system of differential equations y' = f(t,y) by stepping from T0 to 
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution 
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters 
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...). 
%   
%   k_d, gain for the spin rate to reach the sliding mode
%   k_p, gain for the attitude to reach the sliding surface (the bigger it is, the faster the convergence but the bigger the chattering)
 
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.   
%
%   Example 
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);  
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1, 
%     and plots the first component of the solution.   
%
global L F LcHIST y_HIST % global variables control torque vector and torque history
L=zeros(3,1); % L=[L1;L2;L3] total control torque per axis
F=zeros(3,1);
%
if ~isnumeric(tspan)
  error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
  error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
  error('Entries of TSPAN are not in order.') 
end  

try
  f0 = feval(odefun,tspan(1),y0,varargin{:});
catch
  msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
  error(msg);  
end  

y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
  error('Inconsistent sizes of Y0 and f(t0,y0).');
end  

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);

Y(:,1) = y0;

% initialize global variables
LcHIST=zeros(3,N);   % history of control torque vectors
DQHIST=zeros(4,N); % <----- c'est quoi ?????????
DWHIST=zeros(3,N);

% specify static parameters
sample_rate=1;   % sec
%
% principal axis inertia tensor
global Jx Jy Jz

J=[Jx,0,0;0,Jy,0;0,0,Jz];

%
%
% initialize parameters 
% tlast=-557;
% these values for kp and ks will have to be adjusted


y_HIST = zeros(11,N);

for i = 2:N
%
ti = tspan(i-1);
hi = h(i-1);
yi = Y(:,i-1);

y_HIST(:,i-1) = yi;
%
% control algorithm goes here.
% only calculate one control
% distinguish between control sample steps and simulation steps
%
w=yi(5:7);
q3=yi(1:3);
q4=yi(4);
qcs=yi(8:11);

if (ti/sample_rate-floor(ti/sample_rate))<0.001 % test to see if this is a control update (1 msec test threshold)  <---- ????
  % if (ti~=tlast),    % execute if this is not the first sample
    % update control
    % calculate new control

    Xi=[qcs(4),-qcs(3),qcs(2);
        qcs(3),qcs(4),-qcs(1);
        -qcs(2),qcs(1),qcs(4);
        -qcs(1),-qcs(2),-qcs(3)
        ];
    dq3=Xi'*[q3;q4];
    dq4=[q3;q4]'*qcs;

    % this control law normally takes the shortest path to achieve the
    % desired attitude
    L = J*(-k_p*sign(dq4)*dq3 - k_d*(1+dot(dq3', dq3))*w);

    tlast=ti;

else  % this is not a control update sample
% use prior calculated torque, no need to update since global variable
% remembers value
end
% pass torques to odefun
% using global variables now, pass as function parameters possibly eventually
LcHIST(:,i-1)=L(1:3);   % store control torque in global history variable for plotting
DQHIST(1:3,i-1)=dq3'; DQHIST(4,i-1)=dq4;
%
% integrator function calls happen below here
% assumes control torque is the same, even though bfield changes at future times
% if sample time is fast compared to change in bfield, should be okay.
% to improve performance, calculate torque at each intermediate step
% this is where you need interpolation step.
%
 F(:,1) = feval(odefun,ti,yi,varargin{:}); % <-- evaluate dy/dt at ti with y=yi
 F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),varargin{:}); % <-- evaluate dy/dt at ti + 0.5(ti - t(i-1)) with y=yi + 0.5(ti - t(i-1))*dy/dt(ti,yi)
 F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),varargin{:});  
 F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),varargin{:});
 Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
end  % i=2:N

y_HIST;
%
Y = Y.';
