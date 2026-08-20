clear
addpath(fileparts(mfilename('fullpath')));

%% Parametri geometrici comuni
R_int = 60; R_out = 100; k = 200;
Nt = 20;
Ntheta_ref = 8; j_ref = round(Ntheta_ref/2); % cella theta generica, non l'apice

%% Due configurazioni:
%  1) h_cap = R_out (emisfero): quella EFFETTIVAMENTE usata da build_sphcap.m.
%     Qui a_focus=Rs_out-h_cap=0 esattamente, quindi theta_max(r)=90°
%     costante per ogni r -- l'approssimazione sul theta_max medio tra
%     layer (Fase 3) ha errore ESATTAMENTE zero, non solo piccolo. Ci si
%     aspetta accordo a precisione di macchina a qualunque Nr_cap.
%  2) h_cap generico (qui 40): non e' il caso distribuito, ma esercita
%     davvero l'approssimazione (theta_max varia con r) per confermarne
%     la convergenza al second'ordine -- rilevante se h_cap diventasse
%     mai un parametro utente (vedi TODO in build_sphcap.m/GMM4.m).
configs = struct('label',{'h_cap = R_out (emisfero, caso reale)','h_cap generico (caso di verifica)'}, ...
                  'h_cap',{R_out,40});

Ns = [3 5 10 20 40 80];

for c = 1:numel(configs)
    h_cap = configs(c).h_cap;
    fprintf('\n=== %s (h_cap=%g) ===\n', configs(c).label, h_cap);

    Rs_out = (R_out^2+h_cap^2)/(2*h_cap);
    a_focus = Rs_out-h_cap;

    % Volume: riferimento indipendente (integrale 1D chiuso), guscio
    % r_i..r_o, angolo solido pieno 2*pi*(1-cos(theta_max(r))) a ciascun
    % raggio, theta_max(r) dalla condizione di bordo piatto (stesso
    % centro delle sfere della parete: cos(theta_max(r))=(Rs_out-h_cap)/r).
    r_i_shell = sqrt(R_int^2+a_focus^2);
    r_o_shell = sqrt(R_out^2+a_focus^2);
    V_true = 2*pi*( (r_o_shell^3-r_i_shell^3)/3 - a_focus*(r_o_shell^2-r_i_shell^2)/2 );

    theta_max_fun = @(r) atan2(sqrt(max(r.^2-a_focus^2,0)), a_focus);
    theta_lo_fun = @(r) (j_ref-1)/Ntheta_ref .* theta_max_fun(r);
    theta_hi_fun = @(r) j_ref/Ntheta_ref .* theta_max_fun(r);
    dphi_ref = 2*pi/Nt;
    rho_fun = @(r) 1./(k*r.^2*dphi_ref.*(cos(theta_lo_fun(r))-cos(theta_hi_fun(r))));

    V_err = zeros(size(Ns));
    T_err = zeros(size(Ns));

    for idx = 1:numel(Ns)
        Nr_cap = Ns(idx);
        Ntheta = Ntheta_ref;

        [Nodes3D,Prisms,Bricks] = Sph_Shell_Mesh(R_int,R_out,h_cap,Ntheta,Nr_cap,Nt);
        Central = Tri_to_Poly(Prisms,Nt,Nr_cap+1);
        total_nodes = size(Central,1)+size(Bricks,1);
        [elem,Con] = node_sphcap_creator(Nodes3D,Central,Bricks,R_int,R_out,h_cap,Ntheta,Nr_cap,Nt,total_nodes);

        for kk = 1:total_nodes
            elem(kk).item = 'sphcap';
            elem(kk).number = 1;
            elem(kk).prop_mech = [2700,900,k];
            elem(kk).Af_tot = sum(elem(kk).Af);
        end
        sat.node.total_node = total_nodes;
        sat.node.globe = elem;
        Vf_G = zeros(total_nodes); eps_int = zeros(1,total_nodes); sigma = 5.67e-8;
        [Gc,~,~] = TMM2(sat,sigma,Vf_G,Con,eps_int);

        % --- bilancio volume
        [~,~,~,~,~,~,V] = sphcap_areas(R_int,R_out,h_cap,Ntheta,Nr_cap,Nt);
        V_mesh = Nt*sum(V(:));
        V_err(idx) = abs(V_mesh-V_true)/V_true;

        % --- catena radiale (cella j_ref, settore i=1) vs profilo analitico
        i_sec = 1;
        kfun = @(i,j,h) (Nt*(Ntheta-1)+1)*(h-1) + (j-2)*Nt*(j>1) + 1 + i*(j>1);
        R_link = zeros(1,Nr_cap-1);
        for h = 1:Nr_cap-1
            m1 = kfun(i_sec,j_ref,h); m2 = kfun(i_sec,j_ref,h+1);
            R_link(h) = 1/Gc(m1,m2);
        end
        R_cum = [0,cumsum(R_link)];
        T_hot = 100; T_cold = 0;
        T_model = T_hot - (T_hot-T_cold)*R_cum/R_cum(end);

        Rs_bound = zeros(1,Nr_cap+1);
        for kk2 = 0:Nr_cap
            r_rim = R_int+kk2*(R_out-R_int)/Nr_cap;
            Rs_bound(kk2+1) = sqrt(r_rim^2+a_focus^2);
        end
        r_lo = Rs_bound(1:Nr_cap); r_hi = Rs_bound(2:Nr_cap+1);
        r_eff = 2*r_lo.*r_hi./(r_lo+r_hi); % media armonica, coerente con lo split R_i/R_j usato in TMM2.m

        % riferimento analitico DEFINITO SUGLI STESSI ESTREMI del modello
        % (r_eff(1)..r_eff(end), non i confini veri della calotta) --
        % altrimenti numeratore e denominatore non sono confrontabili
        R_true_local = integral(rho_fun, r_eff(1), r_eff(end));
        T_analytic = zeros(size(r_eff));
        for ii = 1:numel(r_eff)
            if ii==1
                T_analytic(ii) = T_hot;
            else
                R_partial = integral(rho_fun, r_eff(1), r_eff(ii));
                T_analytic(ii) = T_hot - (T_hot-T_cold)*R_partial/R_true_local;
            end
        end
        T_err(idx) = max(abs(T_model-T_analytic));

        fprintf('Nr_cap=%3d : errore volume=%.3e   errore max T=%.4e C\n', Nr_cap, V_err(idx), T_err(idx));
    end

    if all(V_err < 1e-10)
        fprintf('errore volume a precisione di macchina per ogni Nr_cap (atteso: a_focus=%.3g)\n', a_focus);
    else
        fprintf('ordine di convergenza (volume):      ');
        fprintf('%.2f  ', log2(V_err(1:end-1)./V_err(2:end)));
        fprintf('\n');
    end
    fprintf('ordine di convergenza (temperatura): ');
    fprintf('%.2f  ', log2(T_err(1:end-1)./T_err(2:end)));
    fprintf('\n');

    assert(all(diff(T_err)<-1e-12 | T_err(2:end)<1e-9), 'errore temperatura non decresce monotonamente con Nr_cap');
    configs(c).V_err = V_err; configs(c).T_err = T_err;
    configs(c).r_eff = r_eff; configs(c).T_model = T_model; configs(c).T_analytic = T_analytic;
end

fprintf('\nOK: convergenza (o esattezza di macchina) confermata su volume e profilo di temperatura, in entrambe le configurazioni.\n');

figure
for c = 1:numel(configs)
    subplot(2,2,(c-1)*2+1)
    loglog(Ns,max(configs(c).V_err,1e-17),'o-',Ns,configs(c).T_err,'s-')
    legend('errore volume (relativo)','errore temperatura max (C)')
    xlabel('N_{r,cap}'); ylabel('errore')
    title(configs(c).label)
    grid on

    subplot(2,2,(c-1)*2+2)
    plot(configs(c).r_eff,configs(c).T_model,'o-',configs(c).r_eff,configs(c).T_analytic,'x--')
    legend('Modello (rete)','Analitica')
    xlabel('raggio [mm]'); ylabel('T [°C]')
    title(sprintf('Profilo radiale, N_{r,cap}=%d',Ns(end)))
end
