function [elem_cap,Con_cap] = build_cyl_cap(R_int,R_out,Nr,Nt,Thickness,rot,Center,z_ref,is_top)
%BUILD_CYL_CAP Full-disc end cap (radius R_out) for a hollow cylinder,
%made of two pieces:
%  - annular layer (R_int..R_out), meshed with the same Nr/Nt resolution
%    as the wall (same Circle_Mesh(R_out,Nr,Nt,R_int) call), so ring
%    boundaries align exactly;
%  - a single lumped node for the inner disc (0..R_int): uniform
%    temperature, capacitance = mass of the whole solid disc. The log
%    formula used for ring-to-ring radial links diverges as r->0, so
%    this region is NOT meshed into further rings; instead the lumped
%    node is linked to the layer's innermost ring (fan-in, all Nt
%    sectors) via a dedicated solid-disc conduction formula in TMM2.m
%    (see the 'lc' node type there).
%Both pieces have axial thickness Thickness (independent of the wall's
%own dz).
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

%% Central lumped node (0..R_int), single isothermal node
% Whole solid disc as one node (not a ring mesh): position at the axis,
% "bottom-face" convention (z=zz(1)) like every other cap element built
% below, later moved to the volumetric center by the shift block. Af(6)
% carries the disc's full exposed base area (pi*R_int^2); Ac is left at
% zero since the radial link to the layer is computed directly in
% TMM2.m, not from Ac/Af areas.
theta = (0:Nt-1)'*(360/Nt);
vertf_local = [R_int*cosd(theta), R_int*sind(theta), repmat(zz(1),Nt,1)];

elem_core = struct( ...
    'ID',[],'ex_in',[],'item',[],'number',[],'ID_item',[], ...
    'node',[0 0 zz(1)]*rot+Center, ...
    'node_diff',[0 0 zz(1)]*rot+Center, ...
    'type','lc', ...
    'face',2, ...
    'vertf',vertf_local*rot+Center, ...
    'Af',[0 0 0 0 0 pi*R_int^2], ...
    'Ac',zeros(1,6), ...
    'V',pi*R_int^2*Thickness, ...
    'prop_mech',[], ...
    'dz_local',Thickness);
Con_core = 0;
n_core = 1;

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

%% Radial stitching: layer's innermost ring (j=1) <-> central lumped node
% Same index formula used internally by node_cyl_creator3 for the layer,
% replicated here since it is not returned by that function. All Nt
% sectors of the layer's bore-facing ring fan into the single lumped
% node (type 'lc'); TMM2.m recognizes the pair via the node type and
% uses the solid-disc conduction formula instead of the log one, which
% diverges as r->0.
k_layer=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;

idx_layer_inner = k_layer(1:Nt,1,1); % layer's bore-facing ring, local indices
m_core = n_layer+1;

for ii=1:1:Nt
    m_layer = idx_layer_inner(ii);

    a_in = elem_cap(m_layer).Af(5); % bore contact area, per sector (external before stitching)
    elem_cap(m_layer).Ac(5) = a_in;
    elem_cap(m_layer).Af(5) = 0;

    Con_cap(m_layer,m_core) = 5; % layer looking inward, toward the core
    Con_cap(m_core,m_layer) = 2; % core looking outward, toward the layer
end

end
