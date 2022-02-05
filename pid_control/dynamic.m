syms q1 q2 q3 q1dot q2dot q3dot q1ddot q2ddot q3ddot;
m_rec = 1e-2;
m_cyl = 5e-2;
c_rec = 0.02*[cos(q1)*sin(q2);sin(q1)*sin(q2);-cos(q2)];
c_cyl = (0.5*q3 + 0.04).*[cos(q1)*sin(q2);sin(q1)*sin(q2);-cos(q2)];
I_cyl = 1e-15*eye(3,3);
I_rec = 1e-15*eye(3,3);
Jv_rec = 0.02*[-sin(q1)*sin(q2) cos(q2)*cos(q1) 0;
              cos(q1)*sin(q2)  cos(q2)*sin(q1) 0;
              0                sin(q2)         0];
r = 0.04 + 0.5*q3;
Jv_cyl = [-r*sin(q1)*sin(q2) r*cos(q2)*cos(q1) cos(q1)*sin(q2);
               r*cos(q1)*sin(q2) r*cos(q2)*sin(q1) sin(q1)*sin(q2);
               0                 r*sin(q2)         -cos(q2)];
Jw_rec = [0 sin(q1) 0; 0 -cos(q1) 0; 1 0 0];
Jw_cyl = [0 sin(q1) 0; 0 -cos(q1) 0; 1 0 0];
R_rec = [cos(q2)*cos(q1) sin(q1) cos(q1)*sin(q2);
         sin(q1)*cos(q2) -cos(q1) sin(q1)*sin(q2);
         sin(q2)          0       -cos(q2)];
R_cyl = R_rec;
%%%%%%%% D-MATRIX COMPUTATION %%%%%%%%%%
D= m_rec*transpose(Jv_rec)*Jv_rec + ...
    transpose(Jw_rec)*R_rec*I_rec*transpose(R_rec)*Jw_rec + ...
    m_cyl*transpose(Jv_cyl)*Jv_cyl + transpose(Jw_cyl)*R_cyl*I_cyl*transpose(R_cyl)*Jw_cyl;
D = vpa(simplify(D),4);
disp('D =');
disp(D)
d11 = D(1,1);
d12 = D(1,2);
d13 = D(1,3);
d21 = D(2,1);
d22 = D(2,2);
d23 = D(2,3);
d31 = D(3,1);
d32 = D(3,2);
d33 = D(3,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% CHRISTOFELL SYMBOL COMPUTATION %%%%
% i = 1, j = 1, k=1
c111 = 0.5*(diff(d11,'q1') + diff(d11,'q1') - diff(d11,'q1'));
% i = 2, j = 1, k=1
c211 = 0.5*(diff(d11,'q2')+diff(d12,'q1') - diff(d21,'q1'));
% i=  3, j=  1, k=1
c311 = 0.5*(diff(d11,'q3')+diff(d13,'q1')- diff(d31,'q1'));
% i= 1, j=2, k=1
c121 = c211;
% i= 2, j=2, k=1
c221 = 0.5*(diff(d12,'q2')+diff(d12,'q2') - diff(d22,'q1'));
% i= 3, j=2, k= 1
c321 = 0.5*(diff(d12,'q3') + diff(d13,'q2') - diff(d32,'q1'));
% i= 1, j=3, k=1
c131 = c311;
% i= 2, j=3, k=1
c231 = c321;
% i= 3, j=3, k=1
c331 = 0.5*(diff(d13,'q3')+diff(d13,'q3')-diff(d33,'q1'));
% i= 1, j=1, k=2
c112 = 0.5*(diff(d21,'q1')+diff(d21,'q1')-diff(d11,'q2'));
% i= 2, j=1, k=2
c212 = 0.5*(diff(d21,'q2')+diff(d22,'q1')-diff(d21,'q2'));
% i= 3, j=1, k=2
c312 = 0.5*(diff(d21,'q3')+diff(d23,'q1')-diff(d31,'q2'));
% i= 1, j=2, k=2
c122 = c212;
% i=2, j=2, k=2
c222 = 0.5*(diff(d22,'q2')+diff(d22,'q2')-diff(d22,'q2'));
% i=3, j=2, k=2
c322 = 0.5*(diff(d22,'q3')+diff(d23,'q2')-diff(d32,'q2'));
% i=1, j=3, k=2
c132 = c312;
% i=2, j=3, k=2
c232 = c322;
% i=3, j=3, k=2
c332 = 0.5*(diff(d23,'q3')+diff(d23,'q3')-diff(d33,'q2'));
% i=1, j=1, k=3
c113 = 0.5*(diff(d31,'q1')+diff(d31,'q1')-diff(d11,'q3'));
% i=2, j=1, k=3
c213 = 0.5*(diff(d31,'q2')+diff(d32,'q1')-diff(d21,'q3'));
% i=3, j=1, k=3
c313 = 0.5*(diff(d31,'q3')+diff(d33,'q1')-diff(d31,'q3'));
% i=1, j=2, k=3
c123 = c213;
% i=2, j=2, k=3
c223 = 0.5*(diff(d32,'q2')+diff(d32,'q2')-diff(d22,'q3'));
% i=3, j=2, k=3
c323 = 0.5*(diff(d32,'q3')+diff(d33,'q2')-diff(d32,'q3'));
% i=1, j=3, k=3
c133 = c313;
% i=2, j=3, k=3
c233 = c323;
% i=3, j=3, k=3
c333 = 0.5*(diff(d33,'q3')+diff(d33,'q3')-diff(d33,'q3'));

C = [c111*q1dot+c211*q2dot+c311*q3dot c121*q1dot+c221*q2dot+c321*q3dot c131*q1dot+c231*q2dot+c331*q3dot; ...
     c112*q1dot+c212*q2dot+c312*q3dot c122*q1dot+c222*q2dot+c322*q3dot c132*q1dot+c232*q2dot+c332*q3dot; ...
     c113*q1dot+c213*q2dot+c313*q3dot c123*q1dot+c223*q2dot+c323*q3dot c133*q1dot+c233*q2dot+c333*q3dot];
C= vpa(simplify(C),4);
disp('C =');
disp(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% GRAVITY TERM COMPUTATION %%%%%%%%%%%%%%

g= [0;9.80665;0];

P = m_rec*g'*c_rec + m_cyl*g'*c_cyl;
g1 = diff(P,'q1');
g2 = diff(P,'q2');
g3 = diff(P,'q3');
g = [g1;g2;g3];
g = vpa(simplify(g),4);
disp('g=');
disp(g);

Tau = D*[q1ddot;q2ddot;q3ddot] + C*[q1dot;q2dot;q3dot] + g;
Tau = vpa(simplify(Tau),4);
disp('Tau=');
disp(Tau);










