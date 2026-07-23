# Changelog — Chiusura del cilindro cavo con tappi meshati (cap_thickness)

Documento di riepilogo delle modifiche per aggiungere ai cilindri cavi la possibilità di
essere chiusi alle due estremità (per modellare un serbatoio invece di un tubo aperto).
Pensato come traccia per la discussione con la relatrice: cosa è stato cambiato, perché,
come è stato verificato, e cosa resta aperto.

Riferimento: lavoro svolto sul branch `TeresaBranch`, in continuazione del lavoro sui
cilindri cavi documentato in `CHANGELOG_cilindro_cavo.md`.

## 1. Problema originale

Un cilindro cavo (`R_int>0`) nel tool è sempre un tubo aperto: il foro centrale non ha mai
nodi, quindi resta aperto su tutta la lunghezza, incluse le due estremità. Per modellare un
serbatoio (parete sottile, chiuso alle due estremità) serve poter "tappare" il foro.

## 2. Prima iterazione (poi sostituita): tappo a nodo singolo

Prima versione, poi abbandonata: un solo nodo isotermo per tappo (area = disco intero
π·R_int², spessore = passo assiale della parete `dz`), collegato radialmente all'anello
più interno della parete con la stessa formula logaritmica già usata per i collegamenti
radiali. Funzionante, ma con una semplificazione forte: nessuna risoluzione radiale
interna al tappo (un solo nodo indipendentemente da `Nr`), e spessore del tappo non
indipendente (sempre uguale a `dz` della parete). Sostituita su richiesta con la versione
meshata descritta sotto.

## 3. Versione finale: tappo meshato, spessore indipendente

Ogni tappo è ora un disco pieno di raggio `R` (raggio esterno della parete), meshato alla
stessa risoluzione `Nr`/`Nt` della parete, con uno spessore fisico proprio (`cap_thickness`,
parametro indipendente da `dz` della parete). Costruito in due pezzi:
- **layer anulare** (da `R_int` a `R`): stessi raggi della parete per costruzione (stessa
  chiamata `Circle_Mesh(R,Nr,Nt,R_int)`), quindi nessuna interpolazione tra mesh diverse;
- **nucleo interno** (da 0 a `R_int`): mini-cilindro pieno (`Circle_Mesh(R_int,Nr,Nt,0)`).

L'altezza `L` inserita dall'utente rappresenta l'altezza totale vera del serbatoio, tappi
inclusi: la parete usa internamente `L - 2·cap_thickness`. Se `L ≤ 2·cap_thickness`, errore
esplicito invece di un fallimento silenzioso.

### File nuovi
- **`Function/build_cyl_cap.m`**: costruisce un tappo completo (layer anulare + nucleo
  interno + loro collegamento radiale). Riusa `Circle_Mesh`/`Mesh2D_to_Mesh3D`/
  `Tri_to_Poly`/`node_cyl_creator3` esattamente come per un cilindro pieno o cavo normale
  — nessuna nuova geometria scritta da zero.
- **`Function/stitch_cyl_wall_and_caps.m`**: combina parete + tappo_bottom + tappo_top in
  un unico `elem`/`Con`, aggiungendo i collegamenti assiali parete↔tappo (stesso i/j, sui
  due anelli affacciati).

### File modificati
- **`Function/node_cyl_creator3.m`**:
  - rimosso il parametro `sealed` e il blocco a nodo singolo (fase 2, superata);
  - aggiunto il campo `dz_local` ad ogni nodo (spessore assiale del *singolo* pezzo di
    mesh costruito in quella chiamata — necessario perché un tappo ha un proprio spessore
    diverso da quello della parete, si veda punto 4);
  - **corretto un bug preesistente**: nel ramo `h==1`, la connessione verso `h+1` veniva
    scritta anche quando quello strato non esiste (`Nz=2`, il caso di un singolo layer,
    usato proprio per il layer anulare del tappo). MATLAB espandeva silenziosamente `Con`
    oltre `total_nodes×total_nodes` invece di dare errore. Innocuo per l'uso storico
    (nessuno chiamava la funzione con `Nz=2`; anche quando succede, le colonne fantasma
    non vengono mai lette da `TMM2.m`), ma bloccante per il tappo meshato. Corretto con una
    guardia `if h<Nz-1`, verificato che per `Nz≥3` (uso normale) è un no-op per costruzione
    logica, non solo empiricamente.
- **`Function/TMM2.m`**: il calcolo della formula logaritmica per i collegamenti radiali
  usava `dz=sat.geom.cyl(cyl_idx).L/(Nz_cyl-1)` — cioè lo spessore assiale **dell'intero
  cilindro** (la parete), sempre lo stesso per qualunque link "same_cyl". Errato per il
  nuovo collegamento radiale *dentro* il tappo, che ha uno spessore fisico diverso
  (`cap_thickness`, non `dz` della parete). Sostituito con `dz` letto dal nuovo campo
  per-nodo `dz_local` di entrambi i nodi del link. Verificato che per ogni nodo esistente
  (parete, cilindro pieno) `dz_local` è bit-per-bit identico al valore che la vecchia
  formula avrebbe calcolato — nessuna differenza per l'uso storico.
- **`Function/GMM4.m`**: legge `cap_thickness` (default 0, retrocompatibile), calcola
  `L_wall=L-2·cap_thickness` con il controllo d'errore, costruisce parete + tappi e li
  unisce con `stitch_cyl_wall_and_caps` quando `R_int>0 && cap_thickness>0`.
- **`Function/node_box_creator2.m`, `node_solid_creator2.m`, `node_face_creator3.m`**:
  aggiunto il campo `dz_local` (vuoto, non usato) alla struct `elem` — necessario solo per
  compatibilità di concatenazione (`GMM4.m` unisce gli `elem` di oggetti diversi in un
  unico array, e MATLAB richiede campi identici tra struct concatenate).

### Bug aggiuntivo corretto durante l'implementazione
Il campo `node_diff`, scritto ma **mai letto da nessuna parte** nel codice, si è rivelato
già malformato in `node_cyl_creator3.m` (manca un `,:` nell'indicizzazione: prende indici
lineari invece di righe, dando un vettore 1×8/1×16 invece di una posizione 1×3). Scoperto
tentando di applicare al tappo la stessa correzione di posizione descritta sotto; non
essendo mai usato altrove, non ha effetto sui risultati — non corretto (fuori scope, non
impatta nulla), solo documentato qui.

### Bug trovato e corretto in fase di test: conduttanza infinita al tappo superiore
Il layer anulare del tappo (mesh a singolo strato, `Nz=2`) eredita da `node_cyl_creator3.m`
la convenzione per cui i nodi di bordo vengono posizionati sulla **superficie esterna**
del layer, non sul baricentro volumetrico. Per il tappo *superiore* questo posizionava il
nodo esattamente alla stessa coordinata z del nodo più esterno della parete — distanza
zero, conduttanza infinita nella formula lineare di `TMM2.m`. Corretto spostando i nodi del
tappo (solo quelli, non quelli della parete) al loro centro volumetrico effettivo
(metà spessore), verificato con test end-to-end che l'errore (NaN/Inf) scompare e le
conduttanze parete↔tappo risultano identiche tra i due lati per simmetria costruttiva.

### Bug trovato e corretto guardando l'anteprima 3D: cilindro fantasma a raggio R_int
Segnalato dalla relatrice guardando il plot: il tappo appariva come se avesse raggio
pari a `R_int`, "infilato dentro" invece che appoggiato sopra. Causa: il nucleo interno
del tappo viene costruito da `build_cyl_cap.m` come un mini-cilindro pieno standalone di
raggio `R_int` — per costruzione, il suo anello più esterno riceve un poligono `vertf` a
due facce (tipo `'s'`: una faccia piatta più una faccia laterale cilindrica), perché per
un cilindro pieno isolato quel bordo *è* davvero la superficie esterna. Nell'assemblaggio
col tappo, invece, quel bordo è un'interfaccia interna (il contatto col layer anulare), non
una superficie esposta — ma il poligono della vecchia faccia laterale (a raggio costante
`R_int`, esteso su tutto lo spessore) restava comunque nella struttura dati e veniva
disegnato da `plot_GMM4.m`, che non fa differenza tra superfici realmente esposte e no.
Le grandezze termiche (`Ac`/`Af`, usate da `TMM2.m`) erano già state corrette in
precedenza per questo bordo — il problema riguardava solo la geometria usata dal plot,
non il calcolo. Corretto in `build_cyl_cap.m`: per questi nodi il tipo viene cambiato da
`'s'` a `'cq'` e il poligono ridotto alla sola faccia piatta (la faccia laterale
fantasma viene scartata). Verificato che dopo la correzione il raggio di tutti i poligoni
disegnabili del cilindro resta entro `[0, R]`, senza più fasce isolate a `R_int`.

## 4. Interfaccia grafica — `TAT2.mlapp`

Nella tabella cilindro (`UITable_5`) la colonna aggiunta in una prima iterazione
(`Closed`, booleano) è stata **rinominata** in `Cap Thickness [mm]` invece di aggiungerne
una nuova: su segnalazione della relatrice, un flag booleano separato sarebbe stato
ridondante rispetto allo spessore stesso (non esiste un caso d'uso sensato per un flag
"chiuso" con spessore zero, né per uno spessore impostato ma ignorato). `cap_thickness`
default `0` → cilindro aperto, comportamento storico invariato.

## 5. Verifica

- **Nessuna regressione**: per `R_int=0` (cilindro pieno) e per `R_int>0` con
  `cap_thickness=0` (cavo aperto, comportamento storico), conteggio nodi e assenza di
  NaN/Inf verificati identici a prima di queste modifiche (132 e 160 nodi rispettivamente,
  stessi valori già osservati in sessioni di validazione precedenti).
- **Caso nuovo** (`R_int=22`, `cap_thickness=5`, stessa mesh Nt=8/Nr=5/Nz=5): 306 nodi
  (160 parete + 73+73 tappi), nessun NaN/Inf, `G_c` simmetrica, rank pieno (306/306).
- Aree dei tappi verificate contro il valore teorico (π·R_int²) tramite
  `center_sort_polygon`/`area_polygon2` (le stesse funzioni usate da `surf_global.m`).
- Test eseguiti sia in isolamento (`build_cyl_cap.m` da solo) sia end-to-end attraverso
  `GMM4.m` + `TMM2.m` con una struttura `sat` completa.

## 6. Limiti noti — non risolti in questo lavoro

1. **Geometria "di faccia" (non di nodo) per plotting/ray-tracing non aggiornata.**
   `cylinder_face.m` costruisce la geometria aggregata delle basi (`sat.geom.cyl(i).face`,
   `sat.geom.globe`) usando ancora l'ingombro della sola parete (`L_wall`), non la vera
   estremità esterna del tappo. Il poligono *per nodo* (`elem.vertf`, quello usato dal
   plot 3D standard `plot_GMM4.m` e da `surf_global.m`) è invece corretto — verificato,
   compreso il fix del cilindro fantasma descritto sopra. Resta da capire se qualche
   altra funzione (es. un ray-tracing a livello di faccia anziché di nodo) usi ancora
   la geometria aggregata non aggiornata; non verificato in questa sessione.
2. **Posizione dei nodi di bordo è un'approssimazione, non il baricentro esatto.** Sia per
   la parete (comportamento preesistente, non toccato) sia per il tappo (dopo la
   correzione del punto 3), i nodi di bordo non sono esattamente al centro del loro
   volume; il collegamento assiale parete↔tappo usa quindi la stessa classe di
   approssimazione già accettata altrove nel tool per i collegamenti assiali tra strati
   della parete — non è stato reso più preciso di così, per coerenza con il resto del
   modello.
3. Come già per `R_int`, questo lavoro riguarda solo la **conduzione**; i limiti noti sul
   modello radiativo del foro documentati in `CHANGELOG_cilindro_cavo.md` restano
   invariati. Il tappo non introduce un modello di contenuto (fluido/propellente): un
   serbatoio con liquido interno richiede un oggetto separato (es. un cilindro concentrico)
   collegato manualmente tramite la tabella "Additional Conduction" — discusso ma non
   implementato in questa sessione.
