function [G0_Irr] = TMM2_only_radiative(sat,sigma,Vf_G,eps_int)

N_nodes=sat.node.total_node;
A_tot_temp = [sat.node.globe.Af_tot];

if size(A_tot_temp,2) > size(A_tot_temp,1)
    A_tot_temp = A_tot_temp';
end
% 10^(-6) multiplication because A_tot is in mm^2
G0_Irr = sigma*10^(-6)*(repmat(A_tot_temp,1,N_nodes).*Vf_G.*repmat(eps_int,1,N_nodes));

end