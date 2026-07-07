function [G_c,G0_Irr,G_hc] = TMM2(sat,sigma,Vf_G,C,eps_int)
nn=sat.node.total_node;

G_c=zeros(nn);                        % Conductor (conduction)
G0_Irr=zeros(nn);                     % Conductor (radiation)
G_hc=zeros(nn);                       % Thermal capacitance

for i=1:1:nn
    G_hc(i,i)=sat.node.globe(i).V*10^(-9)*sat.node.globe(i).prop_mech(2)*...
        sat.node.globe(i).prop_mech(1);
    for j=i:1:nn
        if ne(C(i,j),0)==1
            vect=sat.node.globe(i).node-sat.node.globe(j).node;
            L=norm(vect)/2;
            a1=C(i,j); 
            a2=C(j,i); 
            A1=sat.node.globe(i).Ac(a1);
            A2=sat.node.globe(j).Ac(a2);
            G_i=sat.node.globe(i).prop_mech(3)*A1*(10^-6)/(L*10^-3);
            G_j=sat.node.globe(j).prop_mech(3)*A2*(10^-6)/(L*10^-3);
            G_c(i,j)=(1/G_i+1/G_j)^-1; % series of conductances
            G_c(j,i)=G_c(i,j);
        end  
        if strcmp(sat.node.globe(i).type,'i')==0 && strcmp(sat.node.globe(j).type,'i')==0 % every node but the internal ones
        % eps1=sat.node.globe(i).prop_opt(2);        
        % eps2=sat.node.globe(j).prop_opt(2);
        eps1 = eps_int(i);
        eps2 = eps_int(j);
        G0_Irr(i,j)=sigma*sat.node.globe(i).Af_tot*10^(-6)*Vf_G(i,j)*eps1;
        G0_Irr(j,i)=sigma*sat.node.globe(j).Af_tot*10^(-6)*Vf_G(j,i)*eps2;
        end
    end
end

end








