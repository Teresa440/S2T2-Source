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
            a1=C(i,j);
            a2=C(j,i);

          
            is_radial = (a1==2 && a2==5) || (a1==5 && a2==2);
            same_cyl = strcmp(sat.node.globe(i).item,'cyl') && strcmp(sat.node.globe(j).item,'cyl') ...
                       && sat.node.globe(i).number==sat.node.globe(j).number;

            is_radial_sph = (a1==3 && a2==6) || (a1==6 && a2==3);
            same_sphcap = strcmp(sat.node.globe(i).item,'sphcap') && strcmp(sat.node.globe(j).item,'sphcap') ...
                          && sat.node.globe(i).number==sat.node.globe(j).number;

            use_log = false;
            if is_radial && same_cyl
                [ratio_i,in0_i] = radial_area_ratio(sat.node.globe(i));
                [ratio_j,in0_j] = radial_area_ratio(sat.node.globe(j));
                use_log = ~in0_i && ~in0_j;
            end

            use_sph = is_radial_sph && same_sphcap;

            if use_log
                cyl_idx=sat.node.globe(i).number;
                Nt_cyl=sat.geom.cyl(cyl_idx).Nt;
                alfa=(360/Nt_cyl)/2;

                dz=0.5*(sat.node.globe(i).dz_local+sat.node.globe(j).dz_local);
                C_geom=2*sind(alfa)*dz*10^-3; % [m]: A(r)=C_geom*r for this cylinder's mesh

                R_i=0.5*log(ratio_i)/(sat.node.globe(i).prop_mech(3)*C_geom);
                R_j=0.5*log(ratio_j)/(sat.node.globe(j).prop_mech(3)*C_geom);
                G_c(i,j)=1/(R_i+R_j);
                G_c(j,i)=G_c(i,j);
            elseif use_sph
                % G_r = k*dOmega*r_i*r_o/(r_o-r_i) (Holman Eq. 2-10, hollow
                % sphere), spezzata in due meta' (nodo i, nodo j) come per
                % il ramo cilindrico -- ciascun nodo usa il fattore
                % geometrico della propria cella, gia' incluso dOmega_medio
                % di layer in Geom_r (sphcap_areas.m), qui riusato via
                % dz_local. Convergenza verificata a parte (integrale 1D
                % indipendente, errore ~1/Nr_cap^2).
                Geom_i=sat.node.globe(i).dz_local*10^-3; % [m]
                Geom_j=sat.node.globe(j).dz_local*10^-3; % [m]

                R_i=0.5/(sat.node.globe(i).prop_mech(3)*Geom_i);
                R_j=0.5/(sat.node.globe(j).prop_mech(3)*Geom_j);
                G_c(i,j)=1/(R_i+R_j);
                G_c(j,i)=G_c(i,j);
            else
                % Cucitura parete-calotta: le due mesh condividono
                % esattamente lo stesso confine per costruzione (bordo
                % piatto, Fase 6) -- 'node' e' posizionato apposta sulla
                % faccia di bordo per entrambi i lati e coinciderebbe
                % esattamente (L=0, G_c->Inf). 'node_diff' (baricentro
                % vero, media di tutti gli 8 vertici) e' gia' calcolato
                % da entrambi i node creator per ogni elemento e non
                % soffre di questo problema.
                is_stitch = (strcmp(sat.node.globe(i).item,'cyl') && strcmp(sat.node.globe(j).item,'sphcap')) ...
                         || (strcmp(sat.node.globe(i).item,'sphcap') && strcmp(sat.node.globe(j).item,'cyl'));
                if is_stitch
                    vect=sat.node.globe(i).node_diff-sat.node.globe(j).node_diff;
                else
                    vect=sat.node.globe(i).node-sat.node.globe(j).node;
                end
                L=norm(vect)/2;
                A1=sat.node.globe(i).Ac(a1);
                A2=sat.node.globe(j).Ac(a2);
                G_i=sat.node.globe(i).prop_mech(3)*A1*(10^-6)/(L*10^-3);
                G_j=sat.node.globe(j).prop_mech(3)*A2*(10^-6)/(L*10^-3);
                G_c(i,j)=(1/G_i+1/G_j)^-1; % series of conductances
                G_c(j,i)=G_c(i,j);
            end
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

function [ratio,is_axis] = radial_area_ratio(node)

a_out=node.Ac(2);
if a_out==0
    a_out=node.Af(2);
end
a_in=node.Ac(5);
if a_in==0
    a_in=node.Af(5);
end
is_axis = (a_in==0); % true only for the fused solid-core axis node (r=0)
if is_axis
    ratio = NaN;
else
    ratio = a_out/a_in;
end
end








