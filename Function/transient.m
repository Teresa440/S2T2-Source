function [T,Tempo]=transient(sat,T_0,dt,r_dc,sim_data,Tempo,G0_Irr,G_c,G_hc,...
    Q_0,casevalue)
sigma=5.67e-8;
T=zeros(sat.node.total_node,(length(r_dc)-1)*sim_data.n_orbit);
T(:,1)=T_0;
switch casevalue
    case 1
        diss=sat.node.Q_diss;
        eps_af_earth_matrix=sat.node.analysis.all.epsAf_earth;
        eps_af_space_matrix=sat.node.analysis.all.epsAf_space;
    case 2
        diss=sat.node.Q_diss_cold;
        eps_af_earth_matrix=sat.node.analysis.all.epsAf_earth_cold;
        eps_af_space_matrix=sat.node.analysis.all.epsAf_space_cold;
end

for k=1:(length(Tempo))-1

    Tempo(k+1,1)=Tempo(k,1)+dt;
    epsAf_earth=diag(eps_af_earth_matrix(:,k));
    epsAf_earth_=diag(eps_af_earth_matrix(:,k+1));
    epsAf_space=diag(eps_af_space_matrix(:,k));
    epsAf_space_=diag(eps_af_space_matrix(:,k+1));
    T_=repmat(T(:,k),1,length(T(:,k)));
    T_3=repmat(T(:,k).^3,1,length(T(:,k)));
    T_1_2=(T_.^2)'.*T_;
    T_try_2=T_3+T_3'+T_1_2'+T_1_2;
    G_Irr=G0_Irr'.*T_try_2;                                                 %Gji*(tj-ti)
    G_Irr=G_Irr-diag(sum(G_Irr,2));
    A=(G_c+G_Irr)-diag(4*sigma*epsAf_space_*T(:,k).^3)...
        -diag(4*sigma*epsAf_earth_*T(:,k).^3);
    B=(G_c+G_Irr);
    D=3*sigma*epsAf_space_*(T(:,k).^4)-sigma*epsAf_space*(T(:,k).^4)+...
        Q_0(k,:)'+2*sat.node.Q_diss+Q_0(k+1,:)'+...
        3*sigma*epsAf_earth_*(T(:,k).^4)-...
        sigma*epsAf_earth*(T(:,k).^4);                            %2x Qdiss (t and t+)                           
    AA=G_hc-1/2*A*dt;
    BB=G_hc*T(:,k)+1/2*B*dt*T(:,k)+1/2*D*dt;
    T(:,k+1)=AA\BB;
    toll=1e-3;
    err=100;
    while err>toll
        T_old=T(:,k+1);
        T_=repmat(T_old,1,length(T_old));
        T_3=repmat(T_old.^3,1,length(T_old));
        T_1_2=(T_.^2)'.*T_;
        T_try_2=T_3+T_3'+T_1_2'+T_1_2;
        G_Irr=G0_Irr'.*T_try_2;
        G_Irr=G_Irr-diag(sum(G_Irr,2));
        A=(G_c+G_Irr)-diag(4*sigma*epsAf_space*T_old.^3)...
            -diag(4*sigma*epsAf_earth_*T_old.^3);
        B=(G_c+G_Irr);
        AA=G_hc-1/2*A*dt;
        D=3*sigma*epsAf_space_*(T_old.^4)-sigma*epsAf_space*(T(:,k).^4)+...
        Q_0(k,:)'+2*diss+Q_0(k+1,:)'+...
        3*sigma*epsAf_earth_*(T_old.^4)-...
        sigma*epsAf_earth*(T(:,k).^4); 
        BB=G_hc*T(:,k)+1/2*B*dt*T(:,k)+1/2*D*dt;
        T(:,k+1)=AA\BB;
        err=norm(T(:,k+1)-T_old);
    end
    %constant=1/2*(Q_n+Q_0(k,:))*dt+T(k,:)*G_hc;
    %pause(0.000000000001);
    %app.RunsimulationGauge.Value=20+80*k/(length(Tempo)-1);
end

end