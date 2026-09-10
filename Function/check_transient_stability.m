function check_transient_stability(T, r_dc, n_orbit, node_idx)
% CHECK_TRANSIENT_STABILITY  Verifica la stabilita' numerica dello schema
% Crank-Nicolson di transient.m in risposta alla forzante a gradini
% introdotta da Qdiss_dyn ad ogni transizione di modo operativo.
%
% Controlli:
%  1) Convergenza orbita-su-orbita: max|T_orbit(i+1)-T_orbit(i)| deve
%     tendere a zero al crescere di i (convergenza al regime periodico)
%  2) Ringing sull'ultima orbita: cambi di segno troppo ravvicinati in
%     diff(T), sintomo di oscillazioni numeriche spurie

n_per_orbit = length(r_dc) - 1;

T_node = T(node_idx, :);
T_orbits = reshape(T_node, n_per_orbit, n_orbit)'; % n_orbit x n_per_orbit

%% 1) Convergenza orbita-su-orbita
diff_orbita = zeros(n_orbit-1, 1);
for i = 1:n_orbit-1
    diff_orbita(i) = max(abs(T_orbits(i+1,:) - T_orbits(i,:)));
end

is_monotonic = all(diff(diff_orbita) <= 0);
if ~is_monotonic
    if diff_orbita(end) < 0.1*diff_orbita(1)
        warning('diff_orbita non e'' strettamente monotona ma il trend e'' chiaramente decrescente (oscillazioni smorzate, verosimilmente innocue).');
    else
        warning('diff_orbita non decresce chiaramente verso zero: possibile instabilita'' o mancata convergenza al regime periodico.');
    end
end

figure;
semilogy(1:n_orbit-1, diff_orbita, 'o-', 'LineWidth', 1.3);
xlabel('Indice orbita i (confronto con orbita i+1)');
ylabel('max |T_{i+1} - T_i|  [C]');
title(sprintf('Convergenza orbita-su-orbita -- nodo %d', node_idx));
grid on

%% 2) Ringing sull'ultima orbita
T_last = T_orbits(end,:);
dT = diff(T_last);
sign_dT = sign(dT);
sign_dT(sign_dT == 0) = 1; % evita ambiguita' sui valori esattamente nulli

W = 5;              % ampiezza finestra [campioni]
max_changes_ok = 2; % oltre questo numero di cambi segno nella finestra, e' sospetto

n_flagged = 0;
worst_count = 0;
n_windows = length(sign_dT) - W + 1;
for k = 1:n_windows
    window = sign_dT(k:k+W-1);
    n_changes = sum(diff(window) ~= 0);
    if n_changes > max_changes_ok
        n_flagged = n_flagged + 1;
        worst_count = max(worst_count, n_changes);
    end
end
ringing_detected = n_flagged > 0;

%% Riepilogo
fprintf('\n=== Verifica stabilita'' numerica -- nodo %d ===\n', node_idx);
fprintf('diff_orbita (ultima coppia di orbite): %.4e C\n', diff_orbita(end));
if is_monotonic
    conv_str = 'monotona';
else
    conv_str = 'non monotona (vedi warning sopra)';
end
fprintf('Convergenza orbita-su-orbita: %s\n', conv_str);
fprintf('Finestre con ringing rilevato: %d su %d (finestra=%d campioni)\n', n_flagged, n_windows, W);
fprintf('Massimo numero di cambi di segno in una finestra: %d\n', worst_count);

if diff_orbita(end) > diff_orbita(1)
    giudizio = 'INSTABILE (la differenza tra orbite cresce nel tempo)';
elseif ringing_detected
    giudizio = 'STABILE, ma con rumore numerico locale rilevato (ringing)';
else
    giudizio = 'STABILE';
end
fprintf('Giudizio complessivo: %s\n\n', giudizio);

end
