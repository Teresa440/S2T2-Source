function [r_dc,period,rasc, decl, rsun, sun_vers,v_dc] = Orbit_analysis(env,orbit,sim_data)
%orbita e tempi orbitali
[r_dc,v_dc,period]=orb_calc(env,orbit,sim_data);
%vettore tempo di simulazione
%Tempo=linspace(0,period*sim_data.n_orbit,(period*sim_data.n_orbit)/sim_data.step)';
%Sun Vector
[rasc, decl, rsun] = sun2 (sim_data.jdat);                                  %component sun vector
sun_vers=rsun/norm(rsun);                                                   %sun versor
%angoli di eclisse e tempo di eclisse
[r_dc]=eclipse(r_dc,env,rsun);
%Plot orbit and eclipse period
% plot_mundi(r_dc,env.Rp,2,1);                                              %risolvere il numero di grafico automatico
% figure(2)
% hold on
end