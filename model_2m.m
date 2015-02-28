clear all;
clc;
t=20;
tspan = linspace(0,t,1000*t); 
initialvalue = [-pi/6 0];
[T,Y1]= ode45(@model_2,tspan,initialvalue);
plot(T,Y1(:,1));
ylabel('angle');
xlabel('time');
hold on
plot(T,Y1(:,2),'r')
plot(T,10*pi/3,'g') %required angular velocity
plot(T,pi/6,'-')     %required amplitutde positions
plot(T,-pi/6,'-')
