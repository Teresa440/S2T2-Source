# Changelog — Correzione della conduttanza al nodo d'asse dei cilindri pieni

Documento di riepilogo della correzione applicata alla formula di conduzione usata per il
nodo d'asse fuso (r=0) di **qualunque cilindro pieno** nel tool — non solo il tappo.
Pensato come traccia per la discussione con la relatrice: cosa era sbagliato, come è stato
scoperto, la correzione, e la verifica.

Riferimenti: branch `fix_asse_cilindro_pieno`, creato da `TeresaBranch` (commit `3ae13d6`,
lo stesso punto di partenza di `tappo_singolo_anello`). Il bug e la sua derivazione sono
stati scoperti durante il lavoro documentato in `CHANGELOG_tappo_nodo_lumped.md` (§5 di
quel documento), come verifica collaterale mentre si validava il nodo lumped del tappo —
qui il fix viene applicato al posto giusto: `TMM2.m`, indipendentemente dal tappo.

## 1. Problema

Ogni cilindro pieno (`R_int=0`) nel tool ha, per ciascuno strato assiale, un **nodo
d'asse fuso**: un unico nodo che rappresenta il disco pieno da r=0 al raggio del primo
anello (`a=R/Nr`), collegato all'anello successivo. La formula logaritmica standard usata
per tutti gli altri collegamenti radiali non è definita per r=0 (`ln(r/0)`), quindi
`TMM2.m` la evitava per questo collegamento specifico, usando invece la formula
**lineare/a differenze centrate** generica (la stessa usata per i collegamenti
circonferenziali e assiali, dove è corretta perché lì l'area di conduzione è costante
lungo il percorso).

**Verificato numericamente** (vedi `CHANGELOG_tappo_nodo_lumped.md` §5) che questa
formula lineare applicata al collegamento asse→anello **sottostima la conduttanza di un
fattore ~3**, e — punto chiave — **l'errore non si riduce infittendo la mesh radiale**
(stesso rapporto 0.335 con `Nr=3` e `Nr=50`): non è un errore di discretizzazione, è un
errore strutturale della formula, che tratta un collegamento a sezione convergente verso
r=0 come se fosse a sezione costante.

Questo bug esiste in **tutta la codebase**, per ogni cilindro pieno modellato dal tool —
il tappo di `TeresaBranch` (mesh ad anelli fino al centro) ne è solo un caso particolare,
non l'unico.

## 2. Formula corretta

Stessa derivazione già fatta e validata per il nodo lumped del tappo (vedi
`CHANGELOG_tappo_nodo_lumped.md` §2): conduttanza esatta di un disco pieno di raggio `a`
(qualunque `a`) con generazione volumetrica uniforme fittizia, temperatura del nodo =
media pesata sul volume (coerente con `C=m·cp`):

```
G_disc = 8·π·k·L        (indipendente dal raggio a del disco)
```

ripartita su `Nt` rami paralleli (uno per settore dell'anello a cui il nodo d'asse si
collega):

```
R_core = Nt / (8·π·k·L)
```

in serie con la resistenza del lato anello (**invariata**, stessa formula logaritmica di
sempre).

## 3. Modifica implementata

**`Function/TMM2.m`**: la logica di scelta della formula per i collegamenti radiali
aveva due rami (log per gusci normali, lineare di ripiego per il nodo d'asse). Aggiunto un
terzo ramo, `is_axis_link`, controllato **prima** del ramo log e distinto da esso:

```matlab
is_axis_link = is_radial && same_cyl && (in0_i || in0_j);
```

(`in0_i`/`in0_j` sono lo stesso flag `is_axis` già calcolato da `radial_area_ratio`,
riusato senza modifiche). Quando vero, calcola `R_core` con la formula del punto 2 (usando
`dz_local` e `prop_mech(3)` del nodo d'asse, già disponibili, nessun nuovo campo
necessario) e `R_ring` con la formula log esistente e invariata per il lato anello, poi
`G_c(i,j)=1/(R_ring+R_core)` — stesso schema in serie usato ovunque nel tool.

Il ramo log esistente (`elseif use_log`) e quello lineare generico (per i collegamenti
circonferenziali/assiali, dove restano corretti) **non sono stati toccati**.

**Portata della modifica**: essendo in `TMM2.m` (non in `build_cyl_cap.m` o in un file
specifico del tappo), il fix si applica automaticamente a **qualunque** nodo d'asse di
**qualunque** cilindro pieno gestito dal tool — non solo al tappo che ha innescato la
scoperta del problema.

## 4. Verifica

Tutti i test eseguiti in MATLAB (sessione interattiva, MCP), sul branch
`fix_asse_cilindro_pieno`.

**Test 1 — formula esatta**: isolato il collegamento asse→anello2 di un cilindro pieno di
riferimento (`R=100mm, k=200 W/m·K`), per `Nt∈{8,24}` e `Nr∈{3,8,20,50}`. In tutti gli 8
casi, la conduttanza calcolata dal codice combacia **esattamente** (rapporto 1.000000,
non un'approssimazione) con `8·π·k·dz` calcolata a mano — prima della correzione il
rapporto era 0.335 costante in tutti i casi. La resistenza del lato anello (isolata
sottraendo `R_ring` — calcolata analiticamente con la stessa formula log invariata — dalla
resistenza totale del link) è stata usata per isolare correttamente il solo contributo
`R_core`, dato che con la nuova formula i due lati del collegamento non sono più simmetrici
come lo erano (per costruzione) con la vecchia formula lineare.

**Test 2 — pipeline reale col vecchio tappo** (`build_cyl_cap.m` di `TeresaBranch`, mesh
ad anelli fino al centro, `R=25, R_int=20, L=80, cap_thickness=5, Nr=5, Nt=8, Nz=5`): 306
nodi (160 parete + 73+73 tappi, invariato), `Gc` finita e simmetrica.

**Test 3 — retrocompatibilità, cilindro cavo senza tappo** (`R_int>0`, nessun nodo
d'asse, quindi `is_axis_link` non scatta mai): 640 nodi, nessun NaN/Inf — percorso
strutturalmente isolato dalla modifica, come atteso.

## 5. Impatto sui risultati esistenti — atteso, non un effetto collaterale

Qualunque risultato prodotto in precedenza per un cilindro pieno (con o senza tappo)
**cambierà leggermente vicino all'asse**, perché prima quella regione usava una
conduttanza sbagliata di un fattore ~3. Questo è l'obiettivo della modifica, non un
effetto collaterale indesiderato: i vecchi risultati in quella zona erano scorretti. Se
esistono risultati "validati" storici che si vogliono confrontare, aspettarsi una
differenza proprio nella temperatura vicino al centro dei cilindri pieni, non altrove
(gli altri collegamenti, log e lineare, restano bit-per-bit invariati).

## 6. Nota di scope — relazione con `tappo_singolo_anello`

Questo fix è stato applicato su un branch separato, **a partire da `TeresaBranch`**, non
da `tappo_singolo_anello`. Il nodo lumped del tappo su `tappo_singolo_anello` non aveva
questo problema fin dall'inizio (non usa un nodo d'asse fuso: l'intera regione 0→R_int è
un solo nodo con la formula corretta derivata da zero per quello scopo, si veda
`CHANGELOG_tappo_nodo_lumped.md`) — questo fix non ha quindi nulla da correggere lì, e
`tappo_singolo_anello` non è stato toccato.
