function view_f_plot(sat,r_dc,id_node,surface_n)

figure
plot(r_dc(:,4),sat.node.globe(id_node).view_factor(surface_n).F_p)
hold on
plot(r_dc(:,4),sat.node.globe(id_node).view_factor(surface_n).F_sun)
hold on
plot(r_dc(:,4),r_dc(:,5),'k')
grid on
legend("Planet","Sun","Eclipse")
title("View Factor of node " + string(id_node) + " (type = " + string(sat.node.globe(id_node).item) + ")")
ylabel("F view")
xlabel("Time [s]")

end