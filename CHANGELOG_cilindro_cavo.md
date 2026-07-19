# Changelog — Gestione della geometria dei cilindri cavi (conduzione radiale)

Documento di riepilogo delle modifiche apportate al tool S2T2 per correggere la gestione
della conduzione radiale nei cilindri cavi. Pensato come traccia per la discussione con
la relatrice: cosa è stato cambiato, perché, come è stato verificato, e cosa resta aperto.

Riferimenti commit (branch `TeresaBranch`):
- `cilindro1` — modifiche di codice (geometria, rete di conduttanze, formula di conduzione, script di validazione)
- `cilindro2` — modifica all'interfaccia grafica (TAT2.mlapp)

## 1. Problema originale

Il tool modella un cilindro come una serie di anelli radiali concentrici. Per un cilindro
**pieno**, l'anello più interno tocca l'asse (raggio 0) ed è gestito come caso speciale
(un cuneo/prisma triangolare che converge al centro, fuso in un unico nodo). Il codice
trattava **ogni** cilindro in questo modo, anche quando in realtà avrebbe dovuto
rappresentare un tubo cavo con un foro: l'anello "più interno" veniva sempre trattato come
se toccasse l'asse, anche quando il vero centro del cilindro era vuoto — con area di
conduzione, volume (capacità termica) e posizione del nodo tutti calcolati in modo errato
per quell'anello.

## 2. Modifiche implementate

### Fase 1 — `Function/Circle_Mesh.m`
Aggiunto un parametro opzionale `R_int` (default 0, retrocompatibile). Se `R_int==0`
comportamento identico a prima. Se `R_int>0`, il nodo centrale singolo viene sostituito da
un anello di nodi al raggio `R_int`; i triangoli che convergevano al centro vengono
eliminati e sostituiti da una banda di elementi quadrilateri (`Quads`) che collega
quest'anello al primo anello vero, con lo stesso verso di avvolgimento delle altre bande.

### Fase 2 — `Function/cylinder_areas.m`
Aggiunto `R_int`. Lo spessore radiale (`a`) e le formule d'area per ogni anello ora
partono da `R_int` invece che da 0 (formula della corona circolare invece dello spicchio
pieno). Aggiunta anche la gestione simmetrica della superficie interna del foro: quando
`R_int>0`, la faccia radiale interna del primo anello viene spostata da area di conduzione
(`Ac`) ad area di superficie esterna radiativa (`Af`) — esattamente come già avveniva per
la faccia esterna dell'ultimo anello.

### Fase 3 — `Function/node_cyl_creator3.m`
Quando `R_int>0`, il primo anello (`j==1`) viene trattato come un anello normale (stessa
logica usata per gli anelli interni), non più come colonna centrale fusa. Aggiornata la
formula di indicizzazione dei nodi e dei Bricks per il caso senza nodo centrale, con
guardia per non generare connessioni verso un "anello 0" inesistente quando `R_int>0`.
Il posizionamento del nodo a metà spessore radiale e le aree/volumi corretti derivano
automaticamente dal riuso della stessa logica geometrica già usata per gli anelli
interni — nessun codice ad hoc necessario per quella parte.

**Correzioni necessarie individuate durante l'implementazione** (file coinvolti a catena):
- `Function/Mesh2D_to_Mesh3D.m`: bug latente — l'assegnazione `Triangles(:,4)=-10e10`
  su una matrice vuota (il caso cavo produce `Triangles=[]`) generava una riga fantasma
  scambiata per un "prisma". Corretto con una guardia `if ~isempty(Triangles)`.
- `Function/GMM4.m`: legge il nuovo campo `sat.geom.cyl(i).R_int` (opzionale, default 0
  se assente — retrocompatibile con satelliti già definiti), lo propaga a `Circle_Mesh` e
  `node_cyl_creator3`, e salta la costruzione della colonna centrale (`Tri_to_Poly`) nel
  caso cavo.
- `Function/cylinder_face.m`: il numero di nodi per strato assiale era calcolato come
  `Nt*Nr+1` (valido solo per cilindro pieno); ora calcolato genericamente come
  `size(Nodes3D,1)/Nz`, corretto in entrambi i casi.

### Fase 4 — `Function/TMM2.m`
Per i collegamenti di conduzione radiale (direzione ±r) tra due gusci dello stesso
cilindro, sostituita la formula lineare (approssimazione valida solo per aree di
interfaccia costanti lungo il percorso — vero per le direzioni circonferenziale e
assiale, non per quella radiale) con la formula esatta a due resistori in serie basata
sul logaritmo:

```
R = 0.5 * ln(r_out / r_in) / (k * C_geom)      per ciascun lato del collegamento
```

Il rapporto `r_out/r_in` di ciascun elemento è ricavato da `Ac(2)/Ac(5)` (con fallback su
`Af` quando quella faccia è una superficie fisica esterna — bordo esterno o foro
interno), senza necessità di nuovi campi. Il nodo centrale fuso (solo cilindri pieni,
raggio interno 0) resta sulla vecchia formula lineare, perché il logaritmo non è
definito per r=0 — identificato automaticamente perché lì sia `Ac` che `Af` in direzione
radiale interna sono entrambi zero.

Verifica numerica: per una mesh radiale rada (Nr=4) su un cilindro cavo, la vecchia
formula lineare si discostava dalla nuova fino al 19-23% sui link vicini al foro; la
correzione era quindi non trascurabile proprio nel caso di interesse (mesh rada su un
foro stretto). Il miglioramento vale anche per gli anelli non centrali dei cilindri
pieni, non solo per i cilindri cavi.

### Script di validazione — `Function/validate_hollow_cylinder_conduction.m` (nuovo file)
Script autonomo che costruisce un cilindro cavo di riferimento (R_int=50mm, R_out=100mm,
L=200mm, alluminio k=200 W/mK), esegue l'intera pipeline
(`Circle_Mesh → Mesh2D_to_Mesh3D → node_cyl_creator3 → TMM2`), riduce la rete di
conduttanze a una resistenza termica equivalente e la confronta con la soluzione
analitica esatta di un guscio cilindrico (`R = ln(r2/r1)/(2·π·k·L)`).

Risultato: **errore relativo 0.127%**, coerente con la convergenza attesa per il livello
di raffinamento della mesh usato nel test.

Nota tecnica emersa costruendo il test: la rete di conduttanze collega solo i baricentri
degli anelli, non le superfici fisiche vere (R_int, R_out) — mancano le due "mezze
resistenze" di bordo, che nel modello sono accoppiate solo via radiazione (`Af`), non via
conduzione. Lo script le aggiunge esplicitamente per il confronto bordo-a-bordo (stessa
formula logaritmica di TMM2.m), documentato nei commenti del file.

### Interfaccia grafica — `TAT2.mlapp`
Aggiunta la colonna `R_int [mm]` (15ª colonna, in fondo per non alterare gli indici delle
colonne esistenti) alla tabella di definizione dei cilindri (`UITable_5`):
- `ColumnName`/`ColumnEditable` aggiornati per includere la nuova colonna.
- Valore di default `0` (cilindro pieno, comportamento invariato) nei due blocchi che
  costruiscono le righe della tabella (inizializzazione e aggiunta di nuovi cilindri).
- Il callback che legge la tabella prima di chiamare `GMM4` ora imposta anche
  `sat.geom.cyl(i).R_int` dal valore in tabella.

## 3. Verifica end-to-end (test dal vivo sulla GUI)

Oltre allo script di validazione analitica, la funzionalità è stata testata lanciando
davvero `TAT2`, inizializzando un cilindro dalla UI (R_int=20mm, R=50mm, L=100mm, Nr=5,
Nt=8, Nz=3) e premendo "Preview and Update GMM":
- `sat.geom.cyl(1).R_int`/`.R` propagati correttamente dalla tabella.
- Mesh costruita con 80 nodi = Nt·Nr·(Nz-1), il conteggio atteso per un cilindro cavo
  (nessun nodo centrale fuso).
- Nessun NaN/Inf nelle posizioni dei nodi; raggi dei baricentri correttamente compresi
  tra R_int e R.

## 4. Limiti noti — NON risolti in questo lavoro

Il lavoro sopra riguarda esclusivamente la **conduzione**. Durante l'analisi sono emersi
due gap nel modello **radiativo** (fattori di vista), preesistenti e indipendenti da
queste modifiche, che restano non risolti per scelta (fuori scope, non richiesti):

1. **La parete interna del foro è assente dal calcolo dei fattori di vista.** Gli
   elementi della parete del foro non hanno un ID di superficie assegnato
   (`cylinder_face.m` genera ID solo per le due basi e per la superficie laterale
   esterna, non per quella interna). La funzione `surf_global.m` scarta esplicitamente
   gli elementi senza ID di superficie prima di costruire l'elenco di superfici per il
   ray tracing — quindi quella parete non emette né riceve mai un raggio, in nessuna
   delle due varianti di Monte Carlo ray tracing presenti nel tool (con e senza
   riflessione: entrambe consumano lo stesso elenco di superfici costruito da
   `surf_global.m`).
2. Le basi (corona circolare) e la superficie laterale esterna sono invece risultate
   **corrette** dopo verifica approfondita: un campo "area aggregata" della faccia,
   inizialmente sospettato di ignorare il foro, si è rivelato non utilizzato in nessun
   punto del calcolo effettivo dei fattori di vista (verificato su tutte e quattro le
   varianti di ray tracing presenti nel codice) — quindi non è un bug.

**Rilevanza pratica**: se l'interno del foro non deve scambiare calore per irraggiamento
con nulla (es. contiene un solido/liquido, o è isolato), questo gap non ha effetto sui
risultati. Da notare inoltre che questo tool non ha alcun modello di **convezione** — un
cilindro cavo con liquido all'interno (es. un serbatoio) non è correttamente modellabile
con la sola funzionalità `R_int`, perché il liquido avrebbe massa termica e scambio
convettivo con la parete, entrambi assenti dal modello attuale.
