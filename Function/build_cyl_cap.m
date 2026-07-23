function [elem_cap,Con_cap] = build_cyl_cap(R_int,R_out,Nr,Nt,Thickness,rot,Center,z_ref,is_top)
%BUILD_CYL_CAP Full-disc end cap (radius R_out) for a hollow cylinder,
%meshed at the same Nr/Nt resolution as the wall, made of two pieces:
%  - annular layer (R_int..R_out), built with the same Circle_Mesh(R_out,
%    Nr,Nt,R_int) call as the wall, so ring boundaries align exactly;
%  - inner core (0..R_int), a solid mini-cylinder.
%Both pieces have axial thickness Thickness (independent of the wall's
%own dz) and are stitched together with a radial link (log-formula path
%in TMM2.m, same pattern as a normal ring-to-ring or ring-to-axis link).
%
%z_ref is the wall-facing boundary of the cap, in the wall's *local*
%(pre-rotation) frame: for the bottom cap the disc spans
%[z_ref-Thickness, z_ref]; for the top cap [z_ref, z_ref+Thickness].
%is_top swaps the axial (3<->6) Ac/Af indices and face id, because a
%single-layer (Nz=2) mesh is always built by node_cyl_creator3 as a
%"bottom-type" (h==1) layer.

if is_top
    zz = [z_ref, z_ref+Thickness];
else
    zz = [z_ref-Thickness, z_ref];
end

%% Annular layer (R_int..R_out), same radial grid as the wall
[Nodes_l,Tri_l,Quad_l] = Circle_Mesh(R_out,Nr,Nt,R_int);
[Nodes3D_l,~,Bricks_l] = Mesh2D_to_Mesh3D(Nodes_l,Tri_l,Quad_l,zz);
total_nodes_l = size(Bricks_l,1);
Nodes3D_l = Nodes3D_l*rot+Center;
[elem_layer,Con_layer] = node_cyl_creator3(Nodes3D_l,[],Bricks_l,R_out,Thickness,Nt,Nr,2,total_nodes_l,R_int);
n_layer = numel(elem_layer);

%% Inner core (0..R_int), solid mini-cylinder
[Nodes_c,Tri_c,Quad_c] = Circle_Mesh(R_int,Nr,Nt,0);
[Nodes3D_c,Prisms_c,Bricks_c] = Mesh2D_to_Mesh3D(Nodes_c,Tri_c,Quad_c,zz);
[Central_c] = Tri_to_Poly(Prisms_c,Nt,2);
total_nodes_c = length(Central_c(:,1))+length(Bricks_c(:,1));
Nodes3D_c = Nodes3D_c*rot+Center;
[elem_core,Con_core] = node_cyl_creator3(Nodes3D_c,Central_c,Bricks_c,R_int,Thickness,Nt,Nr,2,total_nodes_c,0);
n_core = numel(elem_core);

%% Top cap: swap axial indices (3<->6), since Nz=2 always yields a
% "bottom-type" (h==1) layer regardless of which physical end it is.
if is_top
    for idx=1:1:n_layer
        elem_layer(idx).Ac([3 6]) = elem_layer(idx).Ac([6 3]);
        elem_layer(idx).Af([3 6]) = elem_layer(idx).Af([6 3]);
        elem_layer(idx).face = 1;
    end
    for idx=1:1:n_core
        elem_core(idx).Ac([3 6]) = elem_core(idx).Ac([6 3]);
        elem_core(idx).Af([3 6]) = elem_core(idx).Af([6 3]);
        elem_core(idx).face = 1;
    end
end

%% Reposition to the cap's own volumetric half-thickness center.
% node_cyl_creator3 places a single-layer (Nz=2) mesh's nodes at zz(1)
% (its "bottom face" convention), i.e. at one edge of the cap, not its
% center. Left as-is, the top cap would land exactly on the wall's own
% h=Nz-1 node (both at z_ref), giving zero distance and Inf conductance
% in TMM2's linear formula. Shifting by half the cap thickness toward
% zz(2) fixes this for both ends, leaving the wall's own node
% convention untouched (consistent with how its existing internal
% axial links already work).
shift_local = [0,0,(zz(2)-zz(1))/2];
shift_global = shift_local*rot;
for idx=1:1:n_layer
    elem_layer(idx).node = elem_layer(idx).node+shift_global;
end
for idx=1:1:n_core
    elem_core(idx).node = elem_core(idx).node+shift_global;
end

%% Combine into one local numbering: layer first, core after
elem_cap = [elem_layer,elem_core];
Con_cap = zeros(n_layer+n_core);
Con_cap(1:n_layer,1:n_layer) = Con_layer;
Con_cap(n_layer+1:end,n_layer+1:end) = Con_core;

%% Radial stitching: layer's innermost ring (j=1) <-> core's outermost ring
% Same index formulas used internally by node_cyl_creator3, replicated
% here since they are not returned by that function.
k_layer=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;                      % R_int>0 branch
k_core=@(i,j,h) (Nt*(Nr-1)+1)*(h-1) + (j-2)*Nt*(j>1) + 1 + i*(j>1); % R_int==0 branch

idx_layer_inner = k_layer(1:Nt,1,1);      % layer's bore-facing ring, local indices
idx_core_outer  = k_core(1:Nt,Nr,1);      % core's outermost ring, local indices (offset by n_layer below)

for ii=1:1:Nt
    m_layer = idx_layer_inner(ii);
    m_core  = idx_core_outer(ii)+n_layer;

    a_in = elem_cap(m_layer).Af(5); % bore contact area, per sector (external before stitching)
    elem_cap(m_layer).Ac(5) = a_in;
    elem_cap(m_layer).Af(5) = 0;

    a_out = elem_cap(m_core).Af(2); % core's outer contact area, per sector (external before stitching)
    elem_cap(m_core).Ac(2) = a_out;
    elem_cap(m_core).Af(2) = 0;

    % The core's outermost ring was built as if it were a standalone
    % solid cylinder of radius R_int, so its 's'-type nodes carry a
    % second polygon (vertf rows 5:8) for what was then a real lateral
    % exterior surface at r=R_int. That boundary is now the internal
    % interface with the annular layer, not an exterior surface -- drop
    % the lateral polygon (keep only the flat wedge face, rows 1:4) so
    % plot_GMM4 doesn't render a phantom small cylinder there.
    if strcmp(elem_cap(m_core).type,'s')
        elem_cap(m_core).type = 'cq';
        elem_cap(m_core).vertf = elem_cap(m_core).vertf(1:4,:);
        elem_cap(m_core).face = elem_cap(m_core).face(1); % drop the lateral face id, keep the top/bottom one
    end

    Con_cap(m_layer,m_core) = 5; % layer looking inward, toward the core
    Con_cap(m_core,m_layer) = 2; % core looking outward, toward the layer
end

end
