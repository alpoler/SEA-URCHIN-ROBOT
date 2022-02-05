% transfer phase trajectory
k=4e-2;
r=4e-2;
t0 = 0;
th = 5;
tf = 10;
foot_zi = -1.5e-2;
foot_zf = 1.5e-2;
foot_yi = 3*sqrt(3)*1e-2;
foot_yh = 4.2e-2;
foot_yf = 3*sqrt(3)*1e-2;
foot_x = 9e-2;

data_mtrx = [1 t0 t0^2 t0^3 t0^4 t0^5;0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;
    0 0 2 6*t0 12*t0^2 20*t0^3;1 tf tf^2 tf^3 tf^4 tf^5;0 1 2*tf 3*tf^2 ...
    4*tf^3 5*tf^4;0 0 2 6*tf 12*tf^2 20*tf^3];
data_mtrxl = [1 t0 t0^2 t0^3 t0^4 t0^5;0 1 2*t0 3*t0^2 4*t0^3 5*t0^4;
    0 0 2 6*t0 12*t0^2 20*t0^3;1 th th^2 th^3 th^4 th^5;0 1 2*th 3*th^2 ...
    4*th^3 5*th^4;0 0 2 6*th 12*th^2 20*th^3];
data_mtrxs = [1 th th^2 th^3 th^4 th^5;0 1 2*th 3*th^2 4*th^3 5*th^4;
    0 0 2 6*th 12*th^2 20*th^3;1 tf tf^2 tf^3 tf^4 tf^5;0 1 2*tf 3*tf^2 ...
    4*tf^3 5*tf^4;0 0 2 6*tf 12*tf^2 20*tf^3];

A1 = data_mtrx \ [foot_zi;0;0;foot_zf;0;0];
Ayl = data_mtrxl \ [foot_yi;0;0;foot_yh;2e-2;0];
Ays = data_mtrxs \ [foot_yh;2e-2;0;foot_yf;0;0];
t = linspace(0,tf,500);
tl = t(1:250);
ts = t(251:500);
z1 = A1(1) + A1(2).*t + A1(3)*(t.^2) + A1(4)*(t.^3) + A1(5)*(t.^4) + A1(6)*(t.^5);
yl = Ayl(1) + Ayl(2).*tl + Ayl(3)*(tl.^2) + Ayl(4)*(tl.^3) + Ayl(5)*(tl.^4) + Ayl(6)*(tl.^5);
ys = Ays(1) + Ays(2).*ts + Ays(3)*(ts.^2) + Ays(4)*(ts.^3) + Ays(5)*(ts.^4) + Ays(6)*(ts.^5);


y = [yl ys];

theta1 = [];
theta2 = [];
d3 = [];


for i=1:500
    [theta_1,theta_2,d_3] = inverse_kinematic(foot_x,y(i),z1(i),r,k);
    theta1 = [theta1 theta_1];
    theta2 = [theta2 theta_2];
    d3 = [d3 d_3];
end




% Support Phase Trajectory

t_sup = linspace(0,3,500);
z_sup = foot_zf - 1e-2.*t_sup;

for i=1:500
    [theta_1,theta_2,d_3] = inverse_kinematic(foot_x,foot_yf,z_sup(i),r,k);

    theta1 = [theta1 theta_1];
    theta2 = [theta2 theta_2];
    d3 = [d3 d_3];
end


t = [t t(end)+t_sup];
theta1v = gradient(theta1,0.02);
theta2v = gradient(theta2,0.02);
d3v = gradient(d3,0.02);
theta1ac = gradient(theta1v,0.02);
theta2ac = gradient(theta2v,0.02);
d3ac = gradient(d3v,0.02);

t = [t t+13 t+26 t+39 t+52];
theta1 = [theta1 theta1 theta1 theta1 theta1];
theta2 = [theta2 theta2 theta2 theta2 theta2];
d3 = [d3 d3 d3 d3 d3];
theta1v = [theta1v theta1v theta1v theta1v theta1v];
theta2v = [theta2v theta2v theta2v theta2v theta2v];
d3v = [d3v d3v d3v d3v d3v];
theta1ac = [theta1ac theta1ac theta1ac theta1ac theta1ac];
theta2ac = [theta2ac theta2ac theta2ac theta2ac theta2ac];
d3ac = [d3ac d3ac d3ac d3ac d3ac];

figure;
subplot(3,3,1);
plot(t,theta1)
title("Joint 1 Angle")
subplot(3,3,2)
plot(t,theta1v)
title("Joint 1 Angular Velocity")
subplot(3,3,3)
plot(t,theta1ac)
title("Joint 1 Angular Accelaration")
subplot(3,3,4)
plot(t,theta2)
title("Joint 2 Angle")
subplot(3,3,5)
plot(t,theta2v)
title("Joint 2 Angular Velocity")
subplot(3,3,6)
plot(t,theta2ac)
title("Joint 2 Angular Acceleration")
subplot(3,3,7)
plot(t,d3)
title("Joint 3 extension")
subplot(3,3,8)
plot(t,d3v)
title("Joint 3 Extension Velocity")
subplot(3,3,9)
plot(t,d3ac)
title("Joint 3 Extension Acceleration")

