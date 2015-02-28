clear all;
[t,y]=ode45(@mod12,[0 10],1)
plot(t,y);
xlabel('time');
ylabel('y');
