clear all;
clc;

options = odeset('RelTol',1e-9,'AbsTol',[1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-5],'MaxStep',0.0005);
t = 25;
tspan = linspace(0,t,1000*t);

k = 0.9;
u = 0.54;

[T,Y1] = ode45(@(t,y) model_1(t,y,u,k),tspan,[0 0 0 0 0 0.00092 400],options);


P = zeros(1000*t,1);
for i = 1:t*1000
    P(i) = 1000*(Y1(i,4)/18+Y1(i,5)/32+Y1(i,6)/28)*8.314*Y1(i,7)/0.001;%Y(i,3)/34
end


% saveas(f,strcat('A_1800_cd_',num2str(u),'_f_',num2str(k),'.jpg'),'jpg');
% %close(f);
% plot(Y(:,1)-Y(:,2));
% 
% % plot(Y(:,2),'r');
% 
% figure;
% plot(Y(:,7));

Y=load('phi 0 10_with pilot.txt','-ascii');
% clc;
% P = load('0_points.txt','-ascii');

title = ['Time'	'Exhaust1'	'Chamber2'	'Pulse'	'Exhaust2'	'Chamber2'	'Exhaust3'	'Chamber3'];
N=length(Y(:,1));
gaussFilter = gausswin(50);
gaussFilter = gaussFilter/sum(gaussFilter);
K = zeros(N,4);
K(:,1) = conv(Y(:,2),gaussFilter,'same');
K(:,2) = conv(Y(:,3),gaussFilter,'same');
K(:,3) = conv(Y(:,4),gaussFilter,'same');
K(:,4) = conv(Y(:,5),gaussFilter,'same');
figure;
plot(Y1(:,7));
f = figure;
plot(P);
% title(strcat('A=1800 cd=',num2str(u),' oxy_fraction=',num2str(k),' \phi =0'));
hold on
plot((K(:,1)+K(:,3)+K(:,4))*100000/3,'r');

% 

