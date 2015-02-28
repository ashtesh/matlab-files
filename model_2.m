function [xdot] = model_2(t,y)
thi = y(1);
z = y(2);
F = 4000;
% F = 5000*(-1*thi/(pi/6));  %newton
if (thi>0) 
    F=(-1)*F;
end
%phi1=0;
r = 0.11;  %m
h = 0.365;
x = 0.11*cos(pi/6);
I = 4;    %kg m2
% phi1 = acos((r-(cos(thi)*x)+(h*sin(thi)))/(sqrt(r^2+x^2+h^2+(2*h*r*sin(thi))-(2*x*r*cos(thi))))) this is dot product;
phi1 = asin(((h*cos(thi))+(x*sin(thi)))/(sqrt(r^2+x^2+h^2+(2*h*r*sin(thi))-(2*x*r*cos(thi)))));%|a cross b|=|a||b||sin(phi)|
% eq =  F*0.11*cos(phi1)/I;
eq = (F*0.11*sin(phi1))/I;
xdot = [z, eq]';
end