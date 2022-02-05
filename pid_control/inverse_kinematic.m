function [ang1,ang2,d3] = inverse_kinematic(fx,fy,fz,r,k)
ang1 = atan2(fy,fx);
ang2 = atan2(sqrt((fy-r)^2 + fx^2),-fz);
d3 = (-fz/cos(ang2)) -r -k;
end
% K corresponds to constant leg length of urchin
% R corresponds to amount of extension of spike after the leg part