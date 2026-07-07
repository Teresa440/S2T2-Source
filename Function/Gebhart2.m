function [Vf_G]=Gebhart2(Vf,eps_for_Ge,sat)

eps=eps_for_Ge;

rif=1-eps;
rif=diag(rif);
eps=diag(eps);
I=eye(sat.node.total_node);
A=(I-Vf*rif);
B=(Vf*eps);
Vf_G=A\B;
end
