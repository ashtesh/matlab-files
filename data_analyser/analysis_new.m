clear; 
 Y=load('comparison_28jan.txt','-ascii');
clc;
 P = load('a_of_0_points_1_28jan.txt','-ascii');

title = ['Time'	'Exhaust1'	'Chamber2'	'Pulse'	'Exhaust2'	'Chamber2'	'Exhaust3'	'Chamber3'];
N=length(Y(:,1));
gaussFilter = gausswin(50);
gaussFilter = gaussFilter/sum(gaussFilter);
K = zeros(N,5);
K(:,1) = conv(Y(:,2),gaussFilter,'same');
K(:,2) = conv(Y(:,3),gaussFilter,'same');
% K(:,3) = conv(Y(:,4),gaussFilter,'same');

k = floor((N)/1000);%k = floor((N-7)/1000);
no = [2 3];    %no = [2 3 4];

inj_pressure = zeros(k,3);
pressure_min = zeros(k,3,2);
inj_prp = zeros(k,3);
inj_67 = zeros(k,3);
p_rate = zeros(k,3);
pmax = zeros(k,3);
p_rate_2 = zeros(k,3);


for a = 1:2;
    for i = 0:89
        inj_pressure(i+1,a) = K(i*1000+1,a);
        [pressure_min(i+1,a,1), pressure_min(i+1,a,2) ]= min(K(i*1000+1:(i+1)*1000,a));
        temp_1 = (inj_pressure(i+1,a)-pressure_min(i+1,a,1))*0.1+pressure_min(i+1,a,1);

%         for j=(pressure_min(i+1,a,2)+i*1000):(i+1)*1000
%             if K(j,a)>inj_pressure(i+1,a);
%                 inj_prp(i+1,a) = j-i*1000-1;
%                 break
%             end
%         end
% 
        for j=(pressure_min(i+1,a,2)+i*1000):(i+1)*1000
            if K(j,a)>temp_1
                inj_67(i+1,a) = j-i*1000-1;
                break
            end
        end
% 
%              
%         pmax (i+1,a) = max(K(i*1000+1:(i+1)*1000,a));
    end
end
% 
% 
% pl = length(P(:,1));
% for i=P(1,1):P(pl-1,1)
%     for a = 1:2
%         xco = (P(i+1,1)*1000+P(i+1,a*2)):(P(i+1,1)*1000+P(i+1,a*2+1));
%         yco = Y((P(i+1,1)*1000+P(i+1,a*2)):(P(i+1,1)*1000+P(i+1,a*2+1)),no(a));
%         slop1 = polyfit(xco,yco',1);
%         p_rate(i+1,a) = slop1(1)*1000;
%         
%     end
% end
%
pl = length(P(:,1));
for i=P(1,1):P(pl-1,1)
    for a = 1:2
        xco = (P(i,1)*1000+P(i+1,a*2)):(P(i,1)*1000+P(i+1,a*2+1));
        yco = Y((P(i,1)*1000+P(i+1,a*2)):(P(i,1)*1000+P(i+1,a*2+1)),no(a));
        slop1 = polyfit(xco,yco',1);
        p_rate(i,a) = slop1(1)*1000;
        
    end
end
% 
% % 
 k_end =89;
 k_start=1;
 %plotprp([k_start:k_end;k_start:k_end;k_start:k_end;k_start:k_end]',[inj_prp(k_start:k_end,1) inj_prp(k_start:k_end,2) inj_prp(k_start:k_end,3) inj_prp(k_start:k_end,4)],'1prp','Chamber',1,'\tau_i_g ,PRP (ms)');
 plotprp([k_start:k_end;k_start:k_end;k_start:k_end]',[inj_67(k_start:k_end,1) inj_67(k_start:k_end,2) inj_67(k_start:k_end,3)],'1.00dp','Chamber',2,'\tau_i_g ,67%\DeltaP (ms)');
% %  plotprp([k_start:k_end;k_start:k_end;k_start:k_end]',[pmax(k_start:k_end,1) pmax(k_start:k_end,2) pmax(k_start:k_end,3)],'1.25pmax','Chamber',3,'P_m_a_x (bar)');
plotprp([k_start:k_end;k_start:k_end;k_start:k_end]',[p_rate(k_start:k_end,1) p_rate(k_start:k_end,2) p_rate(k_start:k_end,3)],'1.0r','Chamber',4,'P_r_a_t_e (bar/sec)');
% % plotprp([k_start:k_end;k_start:k_end;k_start:k_end]',[p_rate_2(k_start:k_end,1) p_rate_2(k_start:k_end,2) p_rate_2(k_start:k_end,3)],'0r','Chamber',4,'P_r_a_t_e 2 (bar/sec)');

% plot(K,'r');
% hold on;
% plot(K(:,1),'r');
% hold on;
% plot(K(:,2),'b');
% plot(Y(:,4),'g');