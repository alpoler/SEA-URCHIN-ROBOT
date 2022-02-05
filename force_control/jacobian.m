function Jmatrix = jacobian(theta1,theta2,d3)
x = (0.04+d3);

Jmatrix = [-x*sin(theta1)*sin(theta2) x*cos(theta2)*cos(theta1) cos(theta1)*sin(theta2); ...
            x*cos(theta1)*sin(theta2) x*cos(theta2)*sin(theta1) sin(theta1)*sin(theta2); ...
            0 x*sin(theta2) -cos(theta2)]';
end
