function [xdot] = model_1(t,x,u,ccc)  
m_in = x(1);
m_out = x(2);
m_h2o2 = x(3);
m_h2o = x(4);
m_o2 = x(5);
m_n2 = x(6);
T = x(7);
%P = x(8);

%% Composition 
    r_h2o2 =0.5 ;
%% Constants
    Ru = 8.314;
%     E =  20400;%J/mole
%     A = 1300; % %*10^6
    E =  40300;%J/mole
    A = 550*260; % %*10^6
    V = 0.001;
    U = 60;
    Ar = 0.04;
%% Density
    d_h2o2 = 1400;
    d_h2o = 1000;
    

    %% Molecular Weight and Mass Fractions
    M_h2o2 = 34;
    M_h2o = 18;
    M_o2 = 32;
    M_n2 = 28;
    
%     m = (m_h2o2 + m_h2o + m_o2 + m_n2);
    m = (m_h2o + m_o2 + m_n2);
%     Y_h2o2 = m_h2o2/m;
    Y_h2o = m_h2o/m;
    Y_o2 = m_o2/m;
    Y_n2 = m_n2/m;
    
%     M = 1/(Y_h2o2/M_h2o2+Y_h2o/M_h2o+Y_o2/M_o2+Y_n2/M_n2);
    M = 1/(Y_h2o/M_h2o+Y_o2/M_o2+Y_n2/M_n2);
    
    R = Ru*1000/M;
    
    
    %% Enthalpy
    h_h2o2 = -5523500; %J/Kg H2O2 Hf0
    h_h2o = -15388000; %J/Kg H2O 
    h_o2 = 0;
    h_n2 = 0;
%     h_in = h_h2o2*r_h2o2+h_h2o*(1-r_h2o2);
%     h_ex = [Y_h2o2, Y_h2o, Y_o2, Y_n2]*[h_h2o2, h_h2o, h_o2, h_n2]';
%% Inlet Exhaust Properties
    C1 = 1;
    C2 = u;
    rin = 0.225;
    A1 = 3.14*10^-6*rin^2;
    r = 1;
    A2 = 3.14*10^-6*r^2;
    P_in = 1100000;
    P_ex = 105000;
%% Heat Capacity J/Kg-K    
    C_h2o2 =  1270+0.618*(T-300);
    C_h2o = 1864+0.62*(T-300);
    C_o2 = 918+0.18*(T-300);
    C_n2 = 1040+0.17*(T-300);
%     C_ex = [C_h2o2, C_h2o, C_o2, C_n2]*[Y_h2o2, Y_h2o, Y_o2, Y_n2]';
%     C = [C_h2o2 C_h2o, C_o2, C_n2]*[m_h2o2 m_h2o m_o2 m_n2]';
    
    C = [C_h2o, C_o2, C_n2]*[m_h2o m_o2 m_n2]';
    C_ex = [ C_h2o, C_o2, C_n2]*[Y_h2o, Y_o2, Y_n2]';
%     mu = Y_h2o*(-4.418944e-06+4.687638e-08*T)+(Y_o2+Y_n2)*(7.879426e-06+4.924946e-08*T);
    %% Equations
 cc = 0.5*square(2*pi()*t,23+15)-0.5*square(2*pi()*t,23);
%    cc = 1;
%     if t<1
%         cc = 0.5+0.5*square(2*pi()*t,ccc);
%     elseif t<2
%         cc = 0.5+0.5*square(2*pi()*t,30);
%     else
%         cc = 0;
%     end
    
    P = 1000*(m_h2o/18+m_o2/32+m_n2/28)*Ru*T/V;%m_h2o2/34
     
    if P_in>P
        m_dot_in = cc*C1*A1*sqrt((d_h2o2*r_h2o2+(1-r_h2o2)*d_h2o)*2*(P_in-P));
    else
        m_dot_in = 0;
    end
    
   
    Pr = P_ex/P;
     C2 = 0.96 - (0.96-0.4)*(P-105000)/325000;
    if Pr<0.528 
          m_dot_out = (A2*C2*P*sqrt((1.4/R)*0.3349))/sqrt(T);
     elseif Pr>1
       m_dot_out = 0;
    else
       m_dot_out = (A2*C2*sqrt((2*1.4/(R*0.4)))*P*(Pr^(1/1.4))*(sqrt(1-(Pr^(0.4/1.4)))))/sqrt(T);
    end
     
    if m_h2o2<0
        m_h2o2=0;
    end

%     if m_h2o2>7e-6
%     rate = A*exp(-E/Ru/T)*7e-6;
%     else
%     rate = A*exp(-E/Ru/T)*m_h2o2;
%     end
        k = 0.9;
        
        rate =0.9*A*exp(-E/Ru/(400))*m_h2o2 + 0.9*A*exp(-E/Ru/600)*m_dot_in*0.002*r_h2o2;

%     rate =0.4*A*exp(-E/Ru/T)*m_h2o2 + A*exp(-E/Ru/T)*m_dot_in*0.0004*r_h2o2;

    m_dot_h2o2 = m_dot_in*r_h2o2 - rate;%- m_dot_out*Y_h2o2;
 
    m_dot_h2o =  rate*(k+M_h2o/M_h2o2) - m_dot_out*Y_h2o;%+m_dot_in*(1-r_h2o2);
 
    m_dot_o2 = rate*(0.5*M_o2/M_h2o2) - m_dot_out*Y_o2;
    
    m_dot_n2 =  - m_dot_out*Y_n2;
    
    dh = 0.5*(2888235 - k*2260000);
    
    T_dot = (rate*dh - m_dot_out*C_ex*(T-298)-0.5*(T-(400-298)*(1-exp(-t/100))-298)+44)/(C-m*R);
    %T_dot = (rate*2888235 - m_dot_out*C_ex*(T-298))/(C-m*R);
        
    %P_dot = 1000*Ru*((m_h2o2/34+m_h2o/18+m_o2/32+m_n2/28)*T_dot+(m_dot_h2o2/34+m_dot_h2o/18+m_dot_o2/32+m_dot_n2/28)*T)/V;
       
    xdot = [m_dot_in, m_dot_out, m_dot_h2o2, m_dot_h2o, m_dot_o2, m_dot_n2, T_dot]';
    
    
end
