%% spin torque efficiency
%[1]2008-Spin-Torque Diode Effect and Its Application-Yoshishige SUZUKI-JPS
%[2]2015-Multiple Reflection Effect on Spin-Transfer Torque-Weiwei Zhu
function f12=torque_eff(mode_select,theta_FL,P1,P2)
switch mode_select

    case 1 %fixed efficiency
        f12=0.8;f21=-0.8;
    case 2 % Slonzwski GMR efficiency

        % to do
    case 3 % Slonzwski TMR efficiency [1]
        f12=P1/(1+P1*P2*cos(theta_FL*pi/180));
    case 4 %multi-reflection :[2]. 2015-Multiple Reflection Effect on Spin-Transfer Torque-Weiwei Zhu
        %% inputs
           %P1,%P2,
        %% outputs
            %f12a%f21a
        %% function defination    
        %[f12a,f21a]=torque_eff(epsilon1_negative,epsilon2_negative,theta_tmp)
        %% start
        %load('torqueeff.mat');
            epsilon1_negative = (1-P1)/2; 
            epsilon2_negative = (1-P2)/2; 
        theta_tmp=theta_FL;
        theta=theta_tmp*pi/180;
        %epsilon1_negative=0.2;
        epsilon1_positive=1-epsilon1_negative;
        %epsilon2_negative=0.2;
        epsilon2_positive=1-epsilon2_negative;
        %theta=pi/4;
        n0_positive=1/sqrt(2);
        n0_negative=1/sqrt(2);

        R1=[1-epsilon1_positive,0;0,1-epsilon1_negative];
        T1=[epsilon1_positive,0;0,epsilon1_negative];
        R2=[1-epsilon2_positive,0;0,1-epsilon2_negative];
        T2=[epsilon2_positive,0;0,epsilon2_negative];

        n0=[n0_positive;n0_negative];

        ot=[(cos(theta/2))^2,(sin(theta/2))^2;...
            (sin(theta/2))^2,(cos(theta/2))^2];%transform operator
        otd=ot';%the conjugate transpose of ot; 'ot dagger'

        I12=(eye(2)-R1*otd*R2*ot)\T1*n0;
        I21=R2*ot*((eye(2)-R1*otd*R2*ot)\T1*n0);

        Iout=otd*T2*ot*((eye(2)-R1*otd*R2*ot)\T1*n0);
        %Iref=T1*otd*R2*ot*((eye(2)-R1*otd*R2*ot)\T1*n0);

        f12=(I12(1)-I12(2))/(Iout(1)+Iout(2));%analytical result 
        f21=(I21(1)-I21(2))/(Iout(1)+Iout(2));%analytical result 
        %f21a=(I21(1)-I21(2))/(Iref(1)+Iref(2));%analytical result 

        %end

end