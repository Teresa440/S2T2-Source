function [E,ni_polar]=inverseproblem(t1p,e,a,mu)
    i=1;
    M=t1p*(mu/a^3)^0.5;
    E(i)=M;
    toll=10e-10;
    err=10e-10;
    while err>=toll
        f(i)=E(i)-e*sin(E(i));
        df_dE(i)=1-e*cos(E(i));                                         %dM/dE
        E(i+1)=E(i)+((M-f(i))/(df_dE(i)));
        err=abs(M-f(i));
        deltaE(i)=abs(E(i+1)-E(i));
        i=i+1;
    end
    
    E_rad=E(end);
    ni_polar=(atan2((1-e^2)^0.5*sin(E_rad),cos(E_rad)-e))*180/pi;
    E=E_rad*180/pi;
    return