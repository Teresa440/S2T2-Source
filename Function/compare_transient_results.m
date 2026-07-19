function compare_transient_results(baseline_path, current_path, do_plot)
%COMPARE_TRANSIENT_RESULTS Confronta due run del transitorio nodo per nodo.
%   compare_transient_results() usa i default 'transient_result_baseline.mat'
%   e 'transient_result.mat' nella cartella corrente.
%   compare_transient_results(baseline_path, current_path) usa percorsi custom.
%   compare_transient_results(baseline_path, current_path, false) disabilita il grafico.

if nargin < 1 || isempty(baseline_path)
    baseline_path = 'transient_result_baseline.mat';
end
if nargin < 2 || isempty(current_path)
    current_path = 'transient_result.mat';
end
if nargin < 3 || isempty(do_plot)
    do_plot = true;
end

B = load(baseline_path, 'sat', 'r_dc');
C = load(current_path, 'sat');

gB = B.sat.node.globe;
gC = C.sat.node.globe;

if numel(gB) ~= numel(gC)
    error('Numero di nodi diverso tra baseline (%d) e run corrente (%d): setup non confrontabile.', ...
        numel(gB), numel(gC));
end

idB = [gB.ID];
idC = [gC.ID];
if ~isequal(idB, idC)
    error('Gli ID dei nodi non coincidono tra baseline e run corrente: setup non confrontabile.');
end

n = numel(gB);
maxAbsErr = zeros(1, n);
maxRelErr = zeros(1, n);

for k = 1:n
    Tb = gB(k).Temperature_t;
    Tc = gC(k).Temperature_t;
    if numel(Tb) ~= numel(Tc)
        error('Nodo %d: lunghezza storia temporale diversa (baseline %d, corrente %d).', ...
            idB(k), numel(Tb), numel(Tc));
    end
    d = abs(Tc(:) - Tb(:));
    maxAbsErr(k) = max(d);
    denom = max(abs(Tb(:)));
    if denom > 0
        maxRelErr(k) = max(d) / denom * 100;
    else
        maxRelErr(k) = 0;
    end
end

[worstAbs, iWorst] = max(maxAbsErr);

fprintf('--- Confronto run transitorio ---\n');
fprintf('Baseline: %s\n', baseline_path);
fprintf('Corrente: %s\n', current_path);
fprintf('Nodi confrontati: %d\n', n);
fprintf('Errore assoluto massimo: %.6g (nodo ID %d)\n', worstAbs, idB(iWorst));
fprintf('Errore relativo massimo: %.6g %%\n', max(maxRelErr));

if worstAbs == 0
    fprintf('Risultato: IDENTICO alla baseline.\n');
else
    fprintf('Risultato: DIVERSO dalla baseline (vedi errori sopra).\n');
end

if ~do_plot
    return
end

Tb = gB(iWorst).Temperature_t(:);
Tc = gC(iWorst).Temperature_t(:);

dt = B.r_dc(2,4) - B.r_dc(1,4);   % same time step the solver uses (Function/transient.m: dt = r_dc(2,4)-r_dc(1,4))
time_h = (0:numel(Tb)-1)' * dt / 3600;   % hours, same axis the software integrates on
diffT = Tc - Tb;

allTb = cat(1, gB.Temperature_t);
rangeGlobal = max(allTb(:)) - min(allTb(:));

figure('Name', 'Transient run comparison', 'Position', [100 100 950 750]);

subplot(2,1,1)
plot(time_h, Tb, 'b-', 'LineWidth', 1.5); hold on
plot(time_h, Tc, 'r--', 'LineWidth', 1.5);
legend('New', 'Original', 'Location', 'best')
xlabel('Time [h]'); ylabel('Temperature [°C]')
title(sprintf('Node ID %d: real scale (total model excursion: %.2f °C)', ...
    idB(iWorst), rangeGlobal))
grid on

subplot(2,1,2)
plot(time_h, diffT, 'k-', 'LineWidth', 1.2);
xlabel('Time [h]'); ylabel('Original - New [°C]')
title(sprintf('Pure difference: max %.4f °C', worstAbs))
grid on
yline(0, '--', 'Color', [0.5 0.5 0.5]);


