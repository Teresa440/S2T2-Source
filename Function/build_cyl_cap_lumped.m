function elem_cap = build_cyl_cap_lumped(R,Nt,cap_thickness,rot,Center,z_ref,is_top)
%BUILD_CYL_CAP_LUMPED Single isothermal end-cap node for a hollow cylinder
%(alternative "lumped" design to build_cyl_cap.m's meshed version): one
%node spanning the *full* cross-section (area = pi*R^2, the whole wall
%radius, not just the bore), volume = pi*R^2*cap_thickness.
%
%This function only builds the node itself (position, exposed area,
%volume, a polygon for plotting). The axial links to the wall's Nr rings
%(one per ring, Con code 7, handled by TMM2.m) are added by the caller,
%since they need the wall's own node indices.
%
%z_ref is the wall-facing boundary of the cap, in the wall's *local*
%(pre-rotation) frame: the cap center sits at z_ref-cap_thickness/2
%(bottom) or z_ref+cap_thickness/2 (top).

if is_top
    z_center_local = z_ref+cap_thickness/2;
    face_id = 1; % top, same convention as the wall's own top face
else
    z_center_local = z_ref-cap_thickness/2;
    face_id = 2; % bottom
end

node_local = [0,0,z_center_local];
node_global = node_local*rot+Center;

alfa = linspace(0,360,Nt+1);
alfa(end) = [];
poly_local = [R*cosd(alfa)',R*sind(alfa)',z_center_local*ones(Nt,1)];
poly_global = poly_local*rot+Center;

elem_cap.ID = [];
elem_cap.ex_in = [];
elem_cap.item = [];
elem_cap.number = [];
elem_cap.ID_item = [];
elem_cap.node = node_global;
elem_cap.node_diff = node_global;
elem_cap.type = 'cap_lumped';
elem_cap.face = face_id;
elem_cap.vertf = poly_global;
elem_cap.Af = zeros(1,6);
if is_top
    elem_cap.Af(3) = pi*R^2;
else
    elem_cap.Af(6) = pi*R^2;
end
elem_cap.Ac = zeros(1,6); % unused for this link type -- see TMM2.m code 7
elem_cap.V = pi*R^2*cap_thickness;
elem_cap.prop_mech = [];
elem_cap.dz_local = cap_thickness;

end
