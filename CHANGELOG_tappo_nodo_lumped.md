# Changelog — Tappo del cilindro cavo: nodo lumped centrale al posto della mesh ad anelli

Documento di riepilogo delle modifiche apportate al tool S2T2 per sostituire, all'interno
del tappo di un cilindro cavo, la mesh ad anelli concentrici fino al centro con un unico
nodo lumped, collegato al resto della mesh con una formula di conduzione derivata
appositamente (non quella logaritmica, che diverge per r→0). Pensato come traccia per la
discussione con la relatrice: cosa è stato cambiato, perché, come è stato derivato e
verificato, e cosa resta aperto.

Riferimenti: lavoro svolto sul branch `tappo_singolo_anello` (creato da `TeresaBranch` al
commit `3ae13d6`, "Tappo meshato per cilindro cavo (cap_thickness)"), commit
`e54e6ee` — "Tappo: nodo lumped centrale al posto della mesh ad anelli fino al centro". In
continuazione diretta del lavoro sul tappo documentato in `CHANGELOG_tappo_cilindro.md` e,
più a monte, su `CHANGELOG_cilindro_cavo.md`.

## 1. Problema originale

Dal commit `3ae13d6`, il nucleo interno del tappo (la regione da r=0 a r=R_int, dentro il
foro della parete cava) era meshato con **anelli concentrici fino al centro**: un
mini-cilindro pieno standalone (`Circle_Mesh(R_int,Nr,Nt,0)` + `node_cyl_creator3`), con
la stessa risoluzione radiale `Nr` del resto del tappo, e un nodo d'asse fuso (r=0) gestito
con la formula lineare di `TMM2.m` (la formula logaritmica, usata per tutti gli altri
collegamenti radiali, non è definita per r=0).

Richiesta della relatrice: sostituire questa regione con un **unico nodo isotermo**
(temperatura uniforme, capacità = massa dell'intero disco pieno interno), mantenendo gli
anelli concentrici solo nella regione anulare R_int→R_out (dove restano coerenti con la
risoluzione della parete). Motivazione: la mesh fine fino al centro non ha senso fisico né
pratico per rappresentare correttamente il tappo — un'unica capacità isoterma è già
un'approssimazione adeguata per quella regione, a patto di collegarla al resto della rete
con una conduttanza corretta, non con la formula lineare che il codice usava già per il
vecchio nodo d'asse (si veda punto 5 più sotto: quella formula si è rivelata sbagliata di
un fattore ~3, indipendentemente dalla mesh).

## 2. Derivazione della formula di conduttanza (nodo lumped ↔ anello)

La formula logaritmica standard per un collegamento radiale fra due gusci cilindrici,

```
R = 0.5 * ln(r_out/r_in) / (k * C_geom)
```

diverge per r_in→0 (già noto e documentato nei changelog precedenti). Per il nodo lumped,
che non è un guscio sottile ma rappresenta l'**intero disco pieno** di raggio `a := R_int`
e spessore assiale `L := cap_thickness`, serve una formula diversa.

### Impostazione

In stazionario puro, senza generazione interna, un disco pieno alimentato solo dal bordo
sarebbe isotermo (nessun flusso netto). Per estrarre una conduttanza equivalente per un
nodo lumped da un solido pieno si assume — come è prassi in letteratura per i "conduttori
a nodo centrale" di cilindri/sfere piene (es. Gilmore, *Spacecraft Thermal Control
Handbook*) — una generazione volumetrica uniforme fittizia `q`, che per bilancio energetico
deve uscire tutta dal bordo r=a. Questo dà un profilo di temperatura parabolico:

```
T(r) = T_bordo + q/(4k)·(a² − r²)
```

con `q = P/(π·a²·L)` (bilancio energetico: tutta la potenza P entra/esce dal bordo).

### Due definizioni possibili di "temperatura del nodo lumped"

**Opzione A — media pesata sul volume** (quella scelta, coerente con la capacità termica
`C=m·cp`, che rappresenta l'energia del nodo come `∫ρ·cp·T dV`, quindi richiede che la
temperatura rappresentativa sia la media volumetrica, non un valore puntuale):

```
T_avg = (2/a²)·∫₀ᵃ T(r)·r dr = T_bordo + q·a²/(8k)
```

Sostituendo `q`:

```
G_A = 8·π·k·L        (indipendente da a)
```

**Opzione B — temperatura puntuale nel centro geometrico** (r=0), scartata:

```
T(0) = T_bordo + q·a²/(4k)  →  G_B = 4·π·k·L   (metà di G_A)
```

L'opzione A è stata confermata esplicitamente dalla relatrice come coerente con la
convenzione già usata nel resto del progetto per le capacità termiche degli altri nodi
(temperatura = media di volume sul volume di controllo, metodo TNM standard).

### Composizione con il lato "anello"

Il collegamento totale anello-più-interno↔nodo-lumped resta in serie, come tutti gli altri
collegamenti radiali del tool:

- **lato anello** (dall'anello più interno del layer verso l'interfaccia r=R_int):
  **invariato**, stessa formula logaritmica già validata, usata così com'è;
- **lato nucleo** (dal nodo lumped verso l'interfaccia r=R_int): nuova formula, `G_A` **
  ripartita su Nt rami paralleli** (uno per ciascun settore angolare dell'anello a cui il
  nodo lumped è collegato):

```
G_core_per_settore = G_A / Nt = 8·π·k·L / Nt
R_core_per_settore = Nt / (8·π·k·L)

G_totale_per_settore = 1 / (R_ring + R_core_per_settore)
```

## 3. Modifiche implementate

### `Function/build_cyl_cap.m`

- **Rimosso** il blocco "Inner core (0..R_int)" che chiamava
  `Circle_Mesh`/`Mesh2D_to_Mesh3D`/`Tri_to_Poly`/`node_cyl_creator3` per costruire un
  mini-cilindro pieno standalone.
- **Sostituito** con la costruzione diretta (nessuna chiamata a mesh-builder) di un
  **singolo elemento** `elem_core`, con un nuovo `type='lc'` ("lumped core") che permette a
  `TMM2.m` di riconoscere il collegamento e instradarlo sulla nuova formula:
  - `node`/`node_diff`: punto sull'asse (r=0), stessa convenzione "faccia inferiore"
    (z=zz(1)) di tutti gli altri elementi del tappo, poi spostato al centro volumetrico dal
    blocco di shift **già esistente e non modificato** (riusato senza toccarlo, valido
    genericamente per un elemento singolo tanto quanto per un array);
  - `face`/`Af`: convenzione "bottom-type" di default (`face=2`, `Af(6)=π·R_int²` — l'intera
    area di base del disco, esposta all'ambiente/irraggiamento), convertita in "top-type"
    (`face=1`, `Af(3)`) dal blocco di swap `is_top` **già esistente, non modificato**;
  - `vertf`: poligono a Nt vertici sul bordo r=R_int (stessa convenzione angolare di
    `cylinder_face.m`), necessario a `surf_global.m` per il ray-tracing Monte Carlo della
    faccia esposta del disco;
  - `Ac`: lasciato a zero — la nuova conduttanza non passa più da rapporti di aree Ac/Af
    (a differenza del collegamento anello-anello), è calcolata direttamente in `TMM2.m`;
  - `V = π·R_int²·Thickness` (volume dell'intero disco pieno, mm³) — la capacità termica
    del nodo (`G_hc`) ne deriva automaticamente, **nessuna modifica necessaria** al calcolo
    generico già presente in `TMM2.m` (`V·prop_mech(2)·prop_mech(1)`);
  - `dz_local = Thickness`, `prop_mech = []` (riempito **automaticamente** dopo, come per
    ogni altro nodo del tappo, da `mech_prop.m`, in base a `item='cyl'`/`number` assegnati
    più tardi da `GMM4.m` — nessuna modifica necessaria lì);
  - **ordine dei campi della struct** replicato esattamente uguale a quello di
    `node_cyl_creator3.m`, condizione necessaria perché MATLAB possa concatenare
    `[elem_layer, elem_core]` e, più a monte, `[elem_wall, elem_cap_bottom, elem_cap_top]`
    in `stitch_cyl_wall_and_caps.m` e infine con gli altri item (`ext`/`sp`/`board`/
    `paral`) in `GMM4.m`.
- **Riscritto** il blocco "Radial stitching" (layer↔nucleo): prima collegava anello-Nr del
  vecchio nucleo standalone all'anello-1 del layer, 1-a-1 su Nt settori (Nt↔Nt); ora
  collega **tutti gli Nt settori dell'anello più interno del layer allo stesso, unico
  nodo lumped** (fan-in Nt↔1). Rimossa anche la logica di ripulitura del poligono
  laterale "fantasma" (`type='s'`→`'cq'`) introdotta nel commit precedente: non più
  necessaria, perché il nodo lumped non è mai stato costruito come cilindro standalone e
  quindi non porta mai quel poligono spurio.

### `Function/TMM2.m`

Aggiunto un terzo ramo alla logica di scelta della formula per i collegamenti radiali
(prima erano solo due: logaritmica per gusci normali, lineare di ripiego per il nodo
d'asse fuso r=0 dei cilindri pieni):

```matlab
is_lc = strcmp(sat.node.globe(i).type,'lc') || strcmp(sat.node.globe(j).type,'lc');
```

Quando `is_radial && same_cyl && is_lc`:
- individua quale dei due nodi è il lumped (`lc`) e quale è l'anello (`ring`);
- calcola `R_core = Nt_cyl/(8·π·k_lc·L_lc)` (formula derivata al punto 2, con
  `L_lc = lc.dz_local` — lo spessore assiale del tappo, lo stesso valore condiviso col
  layer per costruzione);
- calcola `R_ring` con la **stessa identica formula logaritmica già esistente e non
  modificata** (`radial_area_ratio` + mezza-resistenza log), usata finora per il lato
  layer di questo stesso collegamento;
- `G_c(i,j) = 1/(R_ring + R_core)`, `G_c(j,i) = G_c(i,j)`.

Il ramo `is_lc` viene controllato **prima** del ramo log esistente (ora `elseif use_log`)
e lo esclude esplicitamente (`use_log` non viene nemmeno calcolato quando `is_lc` è vero),
per evitare di chiamare `radial_area_ratio` sul nodo lumped — quella funzione assume che
un'area interna nulla significhi "nodo d'asse dei cilindri pieni" (`is_axis`) e andrebbe in
conflitto logico col nuovo caso. Evitato anche il rischio silenzioso opposto: senza questo
ramo dedicato, il nodo lumped (che ha `Ac=zeros(1,6)`) sarebbe caduto nel ramo lineare di
default, con area nulla → conduttanza **zero silenziosa** (non un errore, un bug
subdolo — nessun NaN/Inf, ma un collegamento termico mancante senza avviso).

Nessuna modifica al ramo log esistente, al ramo lineare di default, né al calcolo di
`G_hc`/`G0_Irr` (entrambi già generici rispetto al tipo di nodo).

## 4. Verifica

Tutti i test eseguiti in MATLAB (sessione interattiva, MCP), sia in isolamento sia end-to-end.

**Test A — struttura del nodo lumped** (`R_int=30, R_out=80, Nr=6, Nt=12,
Thickness=4mm, k=200 W/m·K`):
- conteggio nodi: `Nr·Nt+1` (73), confermato;
- `type` dell'ultimo nodo = `'lc'`, confermato;
- `V` calcolato = 11309.733553 mm³, atteso `π·R_int²·Thickness` = 11309.733553 mm³ —
  diff 0;
- `Af(6)` calcolato = 2827.433388 mm², atteso `π·R_int²` = 2827.433388 mm² — match esatto.

**Test B — rete completa via `TMM2.m`**: `Gc` simmetrica (asimmetria massima = 0),
finita ovunque, nessun valore negativo.

**Test C — confronto analitico diretto della conduttanza layer↔nodo lumped**: per
ciascuno dei 12 settori, `Gc(ring1_settore, lc)` confrontata con
`1/(R_ring_analitico + R_core_analitico)` calcolato a mano dagli stessi parametri
geometrici — **errore relativo 0** su tutti i settori (`G = 1.12007962 W/K`).

**Test D — bilancio energetico con dissipazione nota**: rete completa risolta in
stazionario (`A·T = Q`, sistema lineare con eliminazione dei nodi a temperatura fissata),
5 W iniettati nel nodo lumped, anello esterno (j=Nr, tutti i settori) a T=0 fissata.
Confronto fra `T_lumped` del modello e la resistenza composita calcolata a mano (catena
di resistenze log dell'anello + resistenza del nucleo, per simmetria un settore porta
Q/Nt): errore relativo **7.5·10⁻¹⁶** (precisione macchina).

**Test E — retrocompatibilità, percorso cilindro pieno (R_int=0)**: rieseguito lo stesso
tipo di costruzione usata da `Analytical_Validation_cylinder.m` (che non passa mai da
`build_cyl_cap.m`), con `R=100mm, L=200mm, Nr=20, Nt=36, Nz=5, k=200`: 2740 nodi, nessun
NaN/Inf. Il percorso è strutturalmente isolato da queste modifiche (`build_cyl_cap.m`
viene chiamata da `GMM4.m` solo quando `R_int>0 && cap_thickness>0`, si veda
`do_caps` in `GMM4.m`), quindi il confronto conferma che nulla è cambiato lì, non solo
che il risultato è ragionevole.

**Test F — pipeline reale end-to-end** (`node_cyl_creator3` parete + due `build_cyl_cap`
+ `stitch_cyl_wall_and_caps`, esattamente come chiamata da `GMM4.m`, con `R=80, R_int=30,
L=200, Thickness=5, Nr=6, Nt=12, Nz=5`): 288 nodi parete + 73 nodi per tappo (72 layer + 1
lumped) × 2 tappi = 434 nodi totali, esattamente 2 nodi `type='lc'` (uno per tappo), `Gc`
finita e simmetrica sull'intera rete, volumi dei due nodi lumped corretti.

## 5. Scoperta collaterale — la vecchia formula lineare per il nodo d'asse è fisicamente
##    sbagliata di un fattore ~3, indipendentemente dalla mesh

Durante la verifica è emersa una domanda: la vecchia gestione di `TeresaBranch` (nodo
d'asse fuso + formula lineare, invece della formula log che diverge) era perlomeno
fisicamente corretta, anche se meno elegante? Verificato numericamente che **no**.

Isolando il collegamento nodo-d'asse↔anello2 di un cilindro pieno (stesso codice
invariato, presente anche su `TeresaBranch`) e confrontando la conduttanza "disco 0→a"
implicita in quella formula lineare con il valore fisicamente corretto (la stessa
`8·π·k·L` derivata al punto 2, applicata al raggio `a` di quel nodo d'asse):

| Nr (n. anelli) | G formula lineare | G corretta | rapporto |
|---:|---:|---:|---:|
| 3 | 8.4258 W/K | 25.1327 W/K | 0.3353 |
| 8 | 8.4258 W/K | 25.1327 W/K | 0.3353 |
| 20 | 8.4258 W/K | 25.1327 W/K | 0.3353 |
| 50 | 8.4258 W/K | 25.1327 W/K | 0.3353 |

(`R=100mm, Nt=24, Nz=3, k=200 W/m·K`, layer centrale `h=2`)

Il rapporto **non dipende da Nr**: infittire la mesh radiale non riduce l'errore. Facendo
variare invece la risoluzione angolare Nt, il rapporto converge esattamente a **1/3**:

| Nt | rapporto |
|---:|---:|
| 8 | 0.3516 |
| 16 | 0.3377 |
| 24 | 0.3353 |
| 48 | 0.3338 |
| 96 | 0.3335 |

**Interpretazione**: la formula lineare tratta il collegamento asse↔anello2 come una
"barra" a sezione di conduzione costante fra due punti (stesso schema usato per i
collegamenti circonferenziali/assiali, dove è corretto). Vicino a r=0 questo è
strutturalmente sbagliato — non un errore di discretizzazione che si riduce raffinando la
mesh (come sarebbe, ad es., l'errore della formula log applicata a un anello sottile, si
veda `CHANGELOG_cilindro_cavo.md` §Fase 4), ma un errore di **forma della formula**,
costante e persistente qualunque sia `Nr`. Il risultato è un bug silenzioso: nessun
NaN/Inf, temperature apparentemente ragionevoli, ma quantitativamente sbagliate di un
fattore 3 nella regione centrale di ogni cilindro pieno del tool — non solo nel tappo.

**Nota di scope**: questo problema esiste indipendentemente dal lavoro di questo
changelog, in tutta la codebase, per ogni cilindro pieno (`R_int=0`) modellato dal tool,
non solo per i tappi. Non è stato corretto qui (fuori scope, il nodo d'asse dei cilindri
pieni non è stato toccato per garantire la retrocompatibilità richiesta, si veda punto 4
Test E) — segnalato alla relatrice come possibile lavoro futuro.

## 6. Verifica aggiuntiva — lettura diretta dell'ultima run GUI

Ispezionato `GMM_complete.mat`/`Table_Cyl_Geom` dell'ultima run salvata dalla GUI
(`TAT2.mlapp`): `Radius=25mm, R_int=20mm, Length=80mm, Closed(cap_thickness)=5mm, Nr=5,
Nt=8, Nz=5`. Confermato che `Length` è la lunghezza **totale** del cilindro tappi inclusi:
`sat.geom.cyl(1).L=80`, `cap_thickness=5`, e `L_wall = L - 2·cap_thickness = 70mm`
(parete effettiva), coerente con la logica di `GMM4.m` (§3 di
`CHANGELOG_tappo_cilindro.md`). La run in questione conteneva già 2 nodi `type='lc'`,
quindi risulta generata con il codice di questo branch, non con `TeresaBranch`.

## 7. Limiti noti — non risolti in questo lavoro

1. **Fattore ~3 della formula lineare per il nodo d'asse dei cilindri pieni** (punto 5),
   presente in tutta la codebase, non solo nel tappo — fuori scope per questo lavoro,
   segnalato per una futura sessione dedicata.
2. Come per il lavoro precedente sul tappo, restano invariati i limiti noti già
   documentati in `CHANGELOG_tappo_cilindro.md` §6 (geometria "di faccia" aggregata non
   aggiornata per `cylinder_face.m`/`sat.geom.globe`, approssimazione della posizione dei
   nodi di bordo, assenza di un modello di contenuto fluido/propellente nel tappo) e in
   `CHANGELOG_cilindro_cavo.md` §4 (parete interna del foro assente dal calcolo dei
   fattori di vista, nessun modello di convezione).
3. Il nodo lumped, per costruzione, non risolve alcun gradiente radiale interno al disco
   (per definizione — è un unico nodo isotermo). Per casi con transitori molto rapidi e
   `R_int` grande rispetto allo spessore di penetrazione termica, questa è
   un'approssimazione nota e accettata su richiesta esplicita della relatrice (si veda
   punto 1), non un difetto dell'implementazione.
