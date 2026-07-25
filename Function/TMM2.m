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

            % Radial (+-r) links between two shells of the same cylinder have
            % a cross-section that varies with radius along the path: the
            % linear finite-difference formula below is only an approximation
            % there (exact for circumferential/axial links, whose area is
            % constant along the path). Use the exact cylindrical (log) form
            % instead, unless one side is the fused axis node (r=0), where
            % ln(r/0) is undefined and the linear/central formula is kept.
            is_radial = (a1==2 && a2==5) || (a1==5 && a2==2);
            same_cyl = strcmp(sat.node.globe(i).item,'cyl') && strcmp(sat.node.globe(j).item,'cyl') ...
                       && sat.node.globe(i).number==sat.node.globe(j).number;
            is_lc = strcmp(sat.node.globe(i).type,'lc') || strcmp(sat.node.globe(j).type,'lc');

            use_log = false;
            if is_radial && same_cyl && ~is_lc
                [ratio_i,in0_i] = radial_area_ratio(sat.node.globe(i));
                [ratio_j,in0_j] = radial_area_ratio(sat.node.globe(j));
                use_log = ~in0_i && ~in0_j;
            end

            if is_radial && same_cyl && is_lc
                % Cap's central lumped node (whole solid disc, r=0..R_int)
                % <-> annular layer's bore-facing ring: the log formula
                % above diverges as r->0, so the lumped side instead uses
                % the exact steady-state conductance of a solid disc with
                % a uniform (fictitious, energy-balancing) volumetric
                % source, relating its volume-averaged temperature to its
                % boundary temperature at r=R_int: G_disc=8*pi*k*L. Split
                % into Nt parallel branches, one per sector of the ring
                % it fans into. The ring side keeps its usual log
                % half-resistance, unchanged, in series with this one.
                if strcmp(sat.node.globe(i).type,'lc')
                    lc=sat.node.globe(i); ring=sat.node.globe(j);
                else
                    lc=sat.node.globe(j); ring=sat.node.globe(i);
                end
                cyl_idx=ring.number;
                Nt_cyl=sat.geom.cyl(cyl_idx).Nt;
                alfa=(360/Nt_cyl)/2;

                L_disc=lc.dz_local*10^-3; % [m]
                R_core=Nt_cyl/(8*pi*lc.prop_mech(3)*L_disc);

                [ratio_ring,~]=radial_area_ratio(ring);
                C_geom=2*sind(alfa)*ring.dz_local*10^-3; % [m]: A(r)=C_geom*r
                R_ring=0.5*log(ratio_ring)/(ring.prop_mech(3)*C_geom);

                G_c(i,j)=1/(R_ring+R_core);
                G_c(j,i)=G_c(i,j);
            elseif use_log
                cyl_idx=sat.node.globe(i).number;
                Nt_cyl=sat.geom.cyl(cyl_idx).Nt;
                alfa=(360/Nt_cyl)/2;
                % dz_local is per-node (not sat.geom.cyl(cyl_idx).L/Nz),
                % since a link may belong to a cap piece with its own
                % axial thickness, different from the wall's dz.
                dz=0.5*(sat.node.globe(i).dz_local+sat.node.globe(j).dz_local);
                C_geom=2*sind(alfa)*dz*10^-3; % [m]: A(r)=C_geom*r for this cylinder's mesh

                R_i=0.5*log(ratio_i)/(sat.node.globe(i).prop_mech(3)*C_geom);
                R_j=0.5*log(ratio_j)/(sat.node.globe(j).prop_mech(3)*C_geom);
                G_c(i,j)=1/(R_i+R_j);
                G_c(j,i)=G_c(i,j);
            else
                vect=sat.node.globe(i).node-sat.node.globe(j).node;
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
% Outward/inward radial face-area ratio of a cylindrical shell element.
% Ac holds the conduction area unless that face is a real external
% boundary (outer lateral surface, or inner bore for a hollow cylinder),
% in which case the same area is stored in Af instead.
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








