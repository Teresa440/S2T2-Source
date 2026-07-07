function [Central] = Tri_to_Poly(Prisms,Nt,Nz)

Central = zeros(Nz-1,Nt*2);
for h=1:(Nz-1)
    bottom_poly = unique(reshape(Prisms(((h-1)*Nt + 1):(h*Nt), 1:2), 1, [])); % bottom polygons
    top_poly = unique(reshape(Prisms(((h-1)*Nt + 1):(h*Nt), 4:5), 1, [])); % top polygons
    temp_central = [bottom_poly, top_poly];
    Central(h,:) = temp_central;
end

end