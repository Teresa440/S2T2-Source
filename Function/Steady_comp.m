function [T_st]=Steady_comp(sat,Q0_st,G_c,G0_Irr,cases)

sigma=5.67e-8;
err=1e6;
toll=1e-5;
T_st=ones(sat.node.total_node,1)*(273.15+20);
it=1;
switch cases
    case 1
        eps_Af_earth=diag(mean(sat.node.analysis.epsAf_earth,2));
        eps_Af_space=diag(mean(sat.node.analysis.epsAf_space,2));
    case 2
        eps_Af_earth=diag(mean(sat.node.analysis.epsAf_earth_cold,2));
        eps_Af_space=diag(mean(sat.node.analysis.epsAf_space_cold,2));
end

while err > toll
    T_ = repmat(T_st,1,length(T_st));
    T_3 = repmat(T_st.^3,1,length(T_st));
    T_1_2 = (T_.^2)'.*T_;
    T_try_2 = T_3+T_3'+T_1_2'+T_1_2;
    G_Irr = G0_Irr'.*T_try_2;
    G_Irr = G_Irr - diag(sum(G_Irr,2)); % (?)
    A = (G_c + G_Irr) - diag(4*sigma*eps_Af_space*T_st.^3) - diag(4*sigma*eps_Af_earth*T_st.^3); % (?)
    B = -diag(Q0_st) - 3*sigma*eps_Af_space*T_st.^4 - 3*sigma*eps_Af_earth*T_st.^4; % (?)
    T_old = T_st;
    T_st = A\B;
    err = abs(T_old - T_st);
    it = it + 1;
    if it > 1000
        disp('max iteration reached')
        break;
    end
end