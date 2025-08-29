function dydt = eulerseqns_attitude(t,y)
% AE 6356 Assignment 6
% Solution by Glenn Lightsey
% Modified by Alexandre MASSET for Final Project
% November 21, 2024
%
%
% define mass properties
%
global L;
global Jx Jy Jz Jtx Jty Jtz

m = 500;
J=[Jx,0,0;0,Jy,0;0,0,Jz];
Jt=[Jtx,0,0;0,Jty,0;0,0,Jtz];
%
% calculate inertia coefficients
c1=(Jy-Jz)/Jx; 
c2=(Jz-Jx)/Jy; 
c3=(Jx-Jy)/Jz; 
%
% calculate inertia coefficients
ct1=(Jty-Jtz)/Jtx; 
ct2=(Jtz-Jtx)/Jty; 
ct3=(Jtx-Jty)/Jtz; 
%


q_c1=y(1); q_c2= y(2); q_c3=y(3); q_c4=y(4);
w_c1=y(5); w_c2=y(6); w_c3=y(7);
w_c=[w_c1;w_c2;w_c3];
q_t1=y(8); q_t2= y(9); q_t3=y(10); q_t4=y(11);
w_t1=y(12); w_t2=y(13); w_t3=y(14);
w_t=[w_t1;w_t2;w_t3];

skew_w_c = cross_skew_matrix(w_c);
skew_w_t = cross_skew_matrix(w_t);

xi_c=[q_c4,-q_c3,q_c2;
      q_c3,q_c4,-q_c1;
      -q_c2,q_c1,q_c4;
      -q_c1,-q_c2,-q_c3
      ];
xi_t=[q_t4,-q_t3,q_t2;
      q_t3,q_t4,-q_t1;
      -q_t2,q_t1,q_t4;
      -q_t1,-q_t2,-q_t3
      ];
qcd=0.5*xi_c*w_c;
qtd=0.5*xi_t*w_t;

%
% if t<38000
% dydt = [    qcd(1)
%             qcd(2)
%             qcd(3)
%             qcd(4)
%             c1*w_c2*w_c3
%             c2*w_c3*w_c1
%             c3*w_c1*w_c2
%             qtd(1)
%             qtd(2)
%             qtd(3)
%             qtd(4)
%             ct1*w_t2*w_t3 
%             ct2*w_t3*w_t1 
%             ct3*w_t1*w_t2 
%             ];
% else
    dydt = [    qcd(1)
            qcd(2)
            qcd(3)
            qcd(4)
            c1*w_c2*w_c3 + L(1)/Jx
            c2*w_c3*w_c1 + L(2)/Jy
            c3*w_c1*w_c2 + L(3)/Jz
            qtd(1)
            qtd(2)
            qtd(3)
            qtd(4)
            ct1*w_t2*w_t3 
            ct2*w_t3*w_t1 
            ct3*w_t1*w_t2 
            ];
% end

end