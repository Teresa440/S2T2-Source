function [sat] = mech_prop(sat)

for i=1:1:sat.node.total_node
  if strcmp(sat.node.globe(i).item,'ex')==1
      for j=1:1:length(sat.node.globe(i).face)
          face=sat.node.globe(i).face(j);
          sat.node.globe(i).prop_mech=[sat.node.globe(i).prop_mech;sat.prop.ext.mech(face,:)];
      end
      sat.node.globe(i).prop_mech=[mean(sat.node.globe(i).prop_mech(:,1)),...
                                   mean(sat.node.globe(i).prop_mech(:,2)),...
                                   mean(sat.node.globe(i).prop_mech(:,3))];
  elseif strcmp(sat.node.globe(i).item,'sol')==1
      num=sat.node.globe(i).number;
      sat.node.globe(i).prop_mech=sat.prop.sp(num).mech;
  elseif strcmp(sat.node.globe(i).item,'board')==1
      num=sat.node.globe(i).number;
      sat.node.globe(i).prop_mech=sat.prop.board(num).mech;
  elseif strcmp(sat.node.globe(i).item,'paral')==1
      num=sat.node.globe(i).number;
      sat.node.globe(i).prop_mech=sat.prop.parall(num).mech;
  elseif strcmp(sat.node.globe(i).item,'cyl')==1
      num=sat.node.globe(i).number;
      sat.node.globe(i).prop_mech=sat.prop.cyl(num).mech;
  end


end
end