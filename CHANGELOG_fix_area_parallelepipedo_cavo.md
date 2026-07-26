# Changelog — Correzione di area esposta (Af_tot) e volume (V) del parallelepipedo cavo

Documento di riepilogo delle due correzioni applicate al parallelepipedo cavo
(`th>0`): area radiativa esposta (§2-5) e volume/massa/capacità termica (§6-9),
trovate e corrette nella stessa sessione di lavoro. Include anche la ricostruzione
di come questa funzionalità è nata, perché non è stata implementata in questa
sessione né dalla relatrice insieme a me, ma era già presente nel codice quando è
stata controllata.

Riferimenti: branch `fix_af_parallelepipedo_cavo`, creato da `TeresaBranch` (commit
`3ae13d6`, stesso punto di partenza degli altri fix di questa serie).

## 1. Come è nato il parallelepipedo cavo (ricostruito dalla storia git, diff completo)

La funzionalità non è stata introdotta in una sessione di lavoro con me: risulta dal
commit `cb1de4c` ("cilindro_senza_tappo", 23 luglio), **committato direttamente da
Teresa Antico** — nessun trailer `Co-Authored-By: Claude` nel messaggio, a differenza
di tutti i commit dove ho effettivamente lavorato io in questa serie di sessioni.
Genitore del commit: `4323cc1` ("cilindro3"). Il diff completo tra i due tocca 5 file
`.m` (91 righe aggiunte, 7 rimosse) più `TAT2.mlapp` e sei file `.mat` rigenerati; **il
parallelepipedo cavo è solo una delle due funzionalità aggiunte in questo commit**, le
altre righe riguardano un prototipo di tappo per i cilindri, poi abbandonato. Dettaglio
completo, file per file:

### 1.1 `Function/GMM4.m` — la parte rilevante per questo fix (parallelepipedo)

```diff
 L=sat.geom.parall(i).sizes;
 N=sat.geom.parall(i).nodes;
 Angles=sat.geom.parall(i).angles;
+th=sat.geom.parall(i).th;
 ...
-[elem,total_nodes,Connect] = node_solid_creator2(Center,L,N,Angles);
+if th>0
+    Th=repmat(th,1,6);
+    [elem,total_nodes,Connect] = node_box_creator2(Center,L,N,Angles,Th);
+else
+    [elem,total_nodes,Connect] = node_solid_creator2(Center,L,N,Angles);
+end
```

`node_box_creator2.m` è la **stessa identica funzione** già usata per meshare la
Structure esterna (`sat.geom.ext`, con il proprio spessore di parete `th`) — quindi
concettualmente corretta come scelta (un parallelepipedo cavo è geometricamente lo
stesso tipo di oggetto della struttura esterna, solo interno anziché esterno). Il
problema (§2 sotto) è che questo riuso non ha portato con sé una correzione che
esisteva già altrove nel codice proprio per compensare una particolarità di quella
stessa funzione.

### 1.2 `Function/node_cyl_creator3.m` — parte NON rilevante per questo fix (prototipo di tappo per cilindri, abbandonato)

Lo stesso commit aggiunge, allo stesso file, un parametro opzionale `sealed`
(`node_cyl_creator3(...,R_int,sealed)`, default `false`) che — quando vero e
`R_int>0` — aggiunge **due nodi extra** (`type='cap'`) a chiusura del foro, uno per
ciascuna estremità assiale:

```matlab
if sealed && R_int>0
    dz = L/(Nz-1);
    a_in = Afb(1,5);
    Vcap = pi*R_int^2*dz;
    Ac_cap = zeros(1,6); Ac_cap(2) = a_in*Nt;
    Af_cap_bottom(6) = pi*R_int^2;  % area esposta, tappo bottom
    Af_cap_top(3) = pi*R_int^2;     % area esposta, tappo top
    ... % un nodo isotermo per estremità, collegato con Con(...)=2/5
        % a tutti gli Nt settori dell'anello più interno della parete
end
```

Questo è esattamente la **"prima iterazione (poi sostituita): tappo a nodo singolo"**
già descritta in `CHANGELOG_tappo_cilindro.md` §2 — un solo nodo isotermo per tappo,
area = disco intero (`π·R_int²`), collegato radialmente con la stessa formula
logaritmica usata per gli altri collegamenti radiali (non visibile qui, la
formula è scelta da `TMM2.m` in base al codice di collegamento `Con(...)=2/5`, non
codificata in `node_cyl_creator3.m` stesso). Abbandonato e sostituito dal tappo
meshato (`build_cyl_cap.m`) nel commit successivo `3ae13d6`. `Function/GMM4.m` legge
questo flag (`sat.geom.cyl(i).sealed`, default `false` se assente) nella stessa parte
di file, righe adiacenti a quelle del parallelepipedo ma logicamente indipendenti.

### 1.3 Modifiche minori nello stesso commit

- **`Function/Analytical_Validation_cylinder.m`**, **`Function/Validazione2.m`**: una
  riga aggiunta in entrambi, `total_nodes = numel(elem);` dopo la chiamata a
  `node_cyl_creator3`, per contare correttamente i 2 nodi extra quando `sealed=true`
  (difensivo: se `sealed=false`, `numel(elem)` resta invariato).
- **`Function/cylinder_areas.m`**: solo spazi bianchi (trailing spaces aggiunti, una
  riga vuota rimossa) — nessuna modifica funzionale.
- **`TAT2.mlapp`**: modificato (binario, non ispezionabile via diff testuale) —
  presumibilmente per esporre `sealed`/`th` nella UI, non verificato in questa sessione.
- **6 file `.mat`**: rigenerati (dati di un run GUI), non codice.

## 2. Il problema

`node_box_creator2.m` restituisce, per ogni nodo della sua mesh, il campo `Af` (area
radiativa esposta, vettore a 6 direzioni) con l'area reale **duplicata**: gli indici
1&4, 2&5, 3&6 sono sempre uguali a coppie (verificato leggendo il codice sorgente e
poi confermato numericamente, si veda §4). Sommando ingenuamente tutti gli `Af` di
un'intera shell si ottiene quindi **il doppio** dell'area esterna fisica reale.

`Function/opt_prop.m` **conosce già** questo problema, ma solo per l'item Structure:

```matlab
if strcmp(sat.node.globe(i).item,'ex')==1
    sat.node.globe(i).Af_tot=sum(sat.node.globe(i).Af)/2; % only internal-facing area
elseif strcmp(sat.node.globe(i).item,'sol')==1 || strcmp(sat.node.globe(i).item,'board')==1
    sat.node.globe(i).Af_tot=sum(sat.node.globe(i).Af)*2; % top and bottom sides
else
    sat.node.globe(i).Af_tot=sum(sat.node.globe(i).Af);   % <- 'paral' cade qui, nessun /2
end
```

Il ramo `else` (dove cade `item=='paral'`) non applica alcuna correzione. Risultato:
per un parallelepipedo **cavo** (`th>0`, quindi mesh via `node_box_creator2.m`),
`Af_tot` risulta il doppio del valore fisico corretto — e questo entra direttamente
nel calcolo radiativo di `TMM2.m`:

```matlab
G0_Irr(i,j)=sigma*sat.node.globe(i).Af_tot*10^(-6)*Vf_G(i,j)*eps1;
```

quindi lo scambio radiativo di qualunque parallelepipedo cavo verso il resto del
satellite risultava sovrastimato del 100%.

**Non riguarda**: i parallelepipedi pieni (`th=0`, usano `node_solid_creator2.m`,
verificato non avere questa duplicazione — pattern di `Af` con valori distinti, non
a coppie) né la Structure esterna (già corretta con `/2`).

## 3. Correzione implementata

**`Function/GMM4.m`**: invece di estendere la condizione in `opt_prop.m` (che
avrebbe erroneamente dimezzato anche i parallelepipedi pieni, dato che quel ramo è
condiviso da entrambi i casi `th>0`/`th==0`), la correzione è applicata **alla
fonte**, subito dopo la chiamata a `node_box_creator2`, solo nel ramo `th>0`:

```matlab
if th>0
    Th=repmat(th,1,6);
    [elem,total_nodes,Connect] = node_box_creator2(Center,L,N,Angles,Th);
    for j=1:1:total_nodes
        elem(j).Af = elem(j).Af/2;
    end
else
    [elem,total_nodes,Connect] = node_solid_creator2(Center,L,N,Angles);
end
```

Così il ramo generico `else: Af_tot=sum(Af)` di `opt_prop.m` (mai modificato) risulta
già corretto per entrambi i casi, senza toccare `opt_prop.m` né `node_box_creator2.m`
— che restano condivisi con la Structure esterna e con qualunque altro uso futuro,
invariati.

## 4. Verifica

Test eseguiti in MATLAB (sessione interattiva, MCP), sul branch
`fix_af_parallelepipedo_cavo`.

**Verifica preliminare del problema** (prima del fix): costruita una scatola cava di
prova (100×80×60mm, spessore 3mm) con `node_box_creator2.m` isolata, sommati tutti gli
`Af`: **75200 mm²**, contro un'area esterna vera di **37600 mm²** (`2·(Lx·Ly+Ly·Lz+Lx·Lz)`)
— rapporto **2.0000** esatto.

**Test 1 — Af_tot dopo il fix** (stessa geometria, pipeline completa
`GMM4→surf_global→opt_prop`): `Af_tot` totale del parallelepipedo cavo = 37600.00 mm²,
identico all'area vera — rapporto 1.000000.

**Test 2 — retrocompatibilità, parallelepipedo pieno** (`th=0`, stessa geometria):
`Af_tot` = 37600.00 mm², rapporto 1.000000 — invariato, percorso (`node_solid_creator2.m`)
non toccato dal fix.

**Test 3 — retrocompatibilità, Structure esterna**: eseguita senza errori, percorso
(`item=='ex'`, già con la propria correzione `/2` preesistente) non toccato dal fix.

## 5. Impatto sui risultati esistenti (fix dell'area)

Qualunque analisi che abbia usato un parallelepipedo **cavo** (`th>0`) avrà uno
scambio radiativo di quell'oggetto dimezzato rispetto a prima — quello vecchio era
sovrastimato del 100%. Trattandosi di una funzionalità introdotta molto di recente
(commit `cb1de4c`, si veda §1), è probabile che non esistano ancora analisi
"validate" storiche che dipendano dal valore precedente.

## 6. Secondo problema, trovato controllando anche il volume — non solo l'area

Dopo aver corretto `Af`, è stato controllato anche `V` (volume, usato per la massa e
la capacità termica `G_hc` in `TMM2.m`) con lo stesso livello di scrupolo — non per
sospetto specifico, ma perché lo stesso file aveva già dato un bug reale, quindi non
c'era motivo di assumere che il resto fosse corretto senza verificarlo.

**Verifica preliminare**: costruita la stessa scatola cava di prova (100×80×60mm,
spessore 3mm), sommato `V` su tutta la mesh: **112800 mm³**, contro un volume vero del
guscio di **104376 mm³** (`Lx·Ly·Lz − (Lx−2th)(Ly−2th)(Lz−2th)`) — rapporto **1.0807**.

**Per capire se fosse un errore di discretizzazione o un bug vero**: rifatto lo stesso
calcolo a risoluzioni di mesh crescenti (`n=4,6,10,20`, cioè da 56 a 2168 nodi). Il
volume calcolato risulta **identico, bit per bit, a ogni risoluzione** (112800.0000
mm³ sempre) — a differenza, per esempio, dell'approssimazione poligonale del cerchio
nel cilindro (verificata in parallelo per confronto: quella sì converge, da rapporto
0.900 con `Nt=8` a 0.9999 con `Nt=256`). Un valore che non cambia affatto con la
risoluzione, mentre quello vero non dipende dalla mesh, non può essere un effetto di
discretizzazione: è un errore sistematico della formula.

**La causa**: `node_box_creator2.m` calcola il volume di ogni nodo come "superficie
locale del nodo × spessore della propria parete", sommando indipendentemente i
contributi di ciascuna delle pareti che quel nodo tocca (fino a 3, negli spigoli).
Confermato algebricamente: la somma totale coincide **esattamente** con
`area_esterna_scatola × spessore` — la formula "ingenua" che ignora il fatto che, agli
spigoli e agli spifferi dove due o tre pareti si incontrano, quel volume di
sovrapposizione viene contato più volte (due volte su uno spigolo, fino a tre in un
vertice).

**Quanto pesa l'errore**: cresce con lo spessore relativo alle dimensioni della
scatola (non è un valore fisso come per `Af`):

| spessore | errore |
|---|---|
| th=10 (spesso, su una scatola 100×80×60) | +30.6% |
| th=5 | +13.9% |
| th=1 | +2.6% |
| th=0.1 (sottile) | +0.3% |

## 7. Correzione implementata (volume)

Correggere esattamente ogni singolo nodo richiederebbe rifare la geometria di
sovrapposizione per ciascuno dei ~27 casi di `node_box_creator2.m` (spigoli, vertici,
facce), un intervento esteso e a rischio di introdurre nuovi errori proprio nella
funzione dove ne è già stato trovato uno. Scelta invece una correzione **esatta sul
totale, non sul singolo nodo**: nota la formula chiusa del volume vero del guscio,

```matlab
V_true_total = L(1)*L(2)*L(3) - (L(1)-Th(2)-Th(4))*(L(2)-Th(1)-Th(3))*(L(3)-Th(5)-Th(6));
```

(facce 1&3 lungo y, 2&4 lungo x, 5&6 lungo z, secondo le normali di
`face_box_creator.m` — verificate numericamente in questa stessa sessione), si
riscala **ogni** `V` calcolato dello stesso fattore costante:

```matlab
V_naive_total = sum([elem.V]);
for j=1:1:total_nodes
    elem(j).V = elem(j).V*(V_true_total/V_naive_total);
end
```

**Limite consapevole di questa scelta**: il totale (massa e capacità termica
complessiva dell'oggetto) torna esatto, ma la distribuzione **tra i singoli nodi**
resta quella proporzionale della formula originaria (i nodi di spigolo restano quelli
con più volume attribuito, semplicemente tutti scalati della stessa percentuale) —
non è una correzione geometricamente esatta nodo per nodo, a differenza del fix
dell'area (§3), che invece è esatto ovunque, nodo per nodo. Per la stragrande
maggioranza degli usi (bilancio termico complessivo, transitorio dell'oggetto nel suo
insieme) il totale corretto è ciò che conta; un'eventuale analisi che dipendesse dalla
distribuzione fine della capacità termica *dentro* un singolo parallelepipedo cavo
(non tipica per questo genere di modello a rete nodale) risentirebbe ancora di una
leggera imprecisione locale agli spigoli.

## 8. Verifica (volume)

**Test 1 — V_tot esatto dopo il fix**, stessa scatola di prova, spessori
th∈{10,5,1,0.1}: rapporto **1.00000000** esatto in tutti e 4 i casi (prima:
1.3056, 1.1394, 1.0260, 1.0026).

**Test 2 — retrocompatibilità, parallelepipedo pieno** (`th=0`): `V` totale =
480000.0000 mm³ = `Lx·Ly·Lz` esatto — percorso (`node_solid_creator2.m`) non
toccato dal fix.

**Test 3 — verifica che il fix del volume non abbia rotto quello dell'area**:
rieseguito il controllo di `Af_tot` (§4) dopo aver aggiunto anche il fix di `V`:
ancora rapporto 1.00000000 esatto — le due correzioni sono indipendenti, non
interferiscono tra loro.

## 9. Impatto sui risultati esistenti (fix del volume)

Qualunque parallelepipedo cavo già usato in un'analisi avrà ora una massa/capacità
termica minore di prima (tanto meno quanto più lo spessore era grande relativo alle
sue dimensioni) — il transitorio termico di quell'oggetto (velocità di
riscaldamento/raffreddamento) sarà più rapido che nei risultati precedenti. Stessa
considerazione del §5: funzionalità recente, probabilmente nessuna analisi storica
"validata" ne dipende.
