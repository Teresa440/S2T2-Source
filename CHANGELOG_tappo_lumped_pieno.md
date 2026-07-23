# Changelog — Tappo lumped a sezione piena (branch `tappo_lumped_pieno`)

Branch separato da `TeresaBranch` (commit `3ae13d6`, "Tappo meshato per cilindro cavo"),
per provare un'alternativa più semplice ed economica al tappo meshato lì implementato.

## 1. Design

Un solo nodo isotermo per tappo, ma a differenza del primo tentativo lumped (abbandonato
su `TeresaBranch`, area = solo il foro π·R_int²), qui il nodo copre l'**intera sezione**:

- Area (esposta all'ambiente) = **π·R²** (R = raggio esterno della parete, non R_int)
- Volume = **π·R²·cap_thickness**
- Collegato **assialmente** (non radialmente) a *ciascuno* degli `Nr` anelli dell'ultimo
  strato della parete (h=1 per il tappo bottom, h=Nz-1 per il top) — `Nr` collegamenti
  paralleli distinti, uno per anello, non un unico collegamento aggregato.

## 2. File nuovi/modificati

- **`Function/build_cyl_cap_lumped.m`** (nuovo): costruisce il singolo nodo tappo — nessuna
  mesh, nessuna chiamata a `Circle_Mesh`/`node_cyl_creator3`. Include anche un poligono
  (cerchio di raggio R, Nt punti) per una visualizzazione 3D sensata.
- **`Function/TMM2.m`**: aggiunto un nuovo codice di collegamento (**7**, entrambi i lati),
  gestito con una formula dedicata invece della formula lineare generica o di quella
  logaritmica radiale:

  ```
  R_link = (dz_lato1/2)/(k_lato1·A) + (dz_lato2/2)/(k_lato2·A)
  ```

  con `A` = area assiale già calcolata per quell'anello (`Ac(3)`/`Ac(6)`, qualunque dei
  due sia popolato) e `dz_lato/2` = metà spessore di ciascun lato, letto dal campo
  per-nodo `dz_local` (stesso meccanismo già introdotto su `TeresaBranch` per il tappo
  meshato — qui riusato, non reinventato). Non si può usare la formula generica esistente
  perché quella assume una **singola** distanza condivisa (`norm(vect)/2`) uguale per
  entrambi i lati, mentre qui i due spessori (`dz` parete e `cap_thickness`) sono in
  generale diversi.
- **`Function/GMM4.m`**: quando `do_caps` è vero, costruisce i due nodi lumped e, per
  ciascun anello `j=1..Nr` della parete, applica lo stesso scambio `Ac(6)↔Af(6)` (bottom)
  o `Ac(3)↔Af(3)` (top) già usato per il tappo meshato — trasforma la faccia da
  "esposta all'ambiente" a "contatto di conduzione col tappo" — e aggiunge il
  collegamento di codice 7 verso il nodo lumped.

## 3. Cosa NON serve più rispetto al tappo meshato

- Nessuna chiamata a `Circle_Mesh`/`Mesh2D_to_Mesh3D`/`Tri_to_Poly` per il tappo.
- Nessuna divisione layer-anulare/nucleo-interno, nessun collegamento radiale interno
  al tappo, nessuna formula logaritmica per il tappo (`build_cyl_cap.m` e
  `stitch_cyl_wall_and_caps.m` di `TeresaBranch` non sono usati in questo branch).
- Nessun problema di posizione dei nodi da correggere (il nodo lumped non usa la
  posizione automatica per calcolare la distanza di conduzione — la distanza è
  esplicita, `dz_local/2` di ciascun lato — quindi il bug di conduttanza infinita
  trovato sull'altro branch qui non può presentarsi per costruzione).

## 4. Verifica

- **Nodi**: 132 (cilindro pieno, invariato), 160 (cavo aperto, invariato), 162 = 160+2
  (cavo con tappo lumped) — nessuna regressione sui casi esistenti.
- **NaN/Inf**: nessuno. **Simmetria** di `G_c`: esatta. **Rank pieno** (170/170 nel test
  con la struttura esterna inclusa).
- **Conservazione**: somma di riga ~0 dopo la forzatura della diagonale (stesso
  procedimento usato dall'app dopo aver chiamato `TMM2.m`).
- **Verifica indipendente della formula**: calcolato a mano `G` per un singolo
  collegamento anello-tappo (anello j=3, cilindro R=25/R_int=22/cap_thickness=5,
  alluminio k=173) e confrontato con il valore prodotto dalla matrice — **identico**
  (0.1533 W/K in entrambi i casi).
- Poligono del tappo verificato: raggio esattamente R su tutto il perimetro (è un
  cerchio, non una mesh poligonale), volume e area esposta verificati contro le formule
  teoriche (`π·R²·cap_thickness`, `π·R²`).

## 5. Confronto qualitativo coi due approcci

|  | Tappo meshato (`TeresaBranch`) | Tappo lumped pieno (questo branch) |
|---|---|---|
| Risoluzione radiale nel tappo | Sì, stessa `Nr` della parete | No, un solo nodo |
| Costo in nodi aggiunti | ~raddoppia il conteggio nodi del cilindro (dipende da `Nz`) | Solo +2 nodi, costo trascurabile |
| Nuova formula in `TMM2.m` | Riusa quella logaritmica esistente (radiale) | Nuova formula dedicata (codice 7) |
| Complessità implementativa | Alta (3 mesh da incollare, più bug trovati e corretti) | Bassa (nessuna mesh, un nodo e una formula) |
| Realismo fisico | Cattura un eventuale gradiente radiale nel tappo | Assume tappo isotermo (nessun gradiente interno) |

Nessuno dei due è "giusto" in assoluto — la scelta dipende da quanto conta, per lo
studio specifico, risolvere un eventuale gradiente di temperatura dentro il tappo
stesso rispetto al costo computazionale.
