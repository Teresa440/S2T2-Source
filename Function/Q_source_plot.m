function Q_source_plot(r_dc, Source_term, sat)
figure
id_n = 122; % 105; 122;
plot(r_dc(:,4),Source_term.Qa(:,id_n))
hold on
plot(r_dc(:,4),Source_term.Qir(:,id_n))
hold on
plot(r_dc(:,4),Source_term.Qs(:,id_n))
grid on
legend("Planet Albedo","Planet IR","Sun")
title("Heat Sources node " + string(id_n) + " (type = " + string(sat.node.globe(id_n).item) + ")")
ylabel("Q [W]")
xlabel("Time [s]")
end