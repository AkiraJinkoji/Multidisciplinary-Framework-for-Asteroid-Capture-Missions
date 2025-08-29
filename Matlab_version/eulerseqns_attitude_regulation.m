function dydt = eulerseqns_attitude_regulation(t,y)
% AE 6356 Assignment 6
% Solution by Glenn Lightsey
% Modified by Alexandre MASSET for Final Project
% November 21, 2024
%
%
% define mass properties
%
global L;
global Jx Jy Jz

J=[Jx,0,0;0,Jy,0;0,0,Jz];

% calculate inertia coefficients
c1=(Jy-Jz)/Jx; 
c2=(Jz-Jx)/Jy; 
c3=(Jx-Jy)/Jz; 

q_c1=y(1); q_c2= y(2); q_c3=y(3); q_c4=y(4);
w_c1=y(5); w_c2=y(6); w_c3=y(7);
w_c=[w_c1;w_c2;w_c3];

%skew_w_c = cross_skew_matrix(w_c);

xi_c=[q_c4,-q_c3,q_c2;
      q_c3,q_c4,-q_c1;
      -q_c2,q_c1,q_c4;
      -q_c1,-q_c2,-q_c3
      ];

qcd=0.5*xi_c*w_c;

dydt = [    qcd(1)
        qcd(2)
        qcd(3)
        qcd(4)
        c1*w_c2*w_c3 + L(1)/Jx
        c2*w_c3*w_c1 + L(2)/Jy
        c3*w_c1*w_c2 + L(3)/Jz
        0
        0
        0
        0
        ];


end