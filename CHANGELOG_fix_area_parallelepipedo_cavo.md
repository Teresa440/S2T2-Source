# Changelog — Correzione dell'area esposta (Af_tot) del parallelepipedo cavo

Documento di riepilogo della correzione applicata all'area radiativa esposta di un
parallelepipedo cavo (`th>0`). Include anche la ricostruzione di come questa
funzionalità è nata, perché non è stata implementata in questa sessione né dalla
relatrice insieme a me, ma era già presente nel codice quando è stata controllata.

Riferimenti: branch `fix_af_parallelepipedo_cavo`, creato da `TeresaBranch` (commit
`3ae13d6`, stesso punto di partenza degli altri fix di questa serie).

## 1. Come è nato il parallelepipedo cavo (ricostruito dalla storia git)

La funzionalità non è stata introdotta in una sessione di lavoro con me: risulta dal
commit `cb1de4c` ("cilindro_senza_tappo", 23 luglio), **committato direttamente da
Teresa Antico** — nessun trailer `Co-Authored-By: Claude` nel messaggio, a differenza
di tutti i commit dove ho effettivamente lavorato io in questa serie di sessioni. Lo
stesso commit introduceva anche un parametro `sealed` per i cilindri (un prototipo
precoce di quello che sarebbe poi diventato `R_int`, poi superato dal lavoro
documentato in `CHANGELOG_cilindro_cavo.md`).

La modifica in `Function/GMM4.m` per il parallelepipedo è minimale e riusa codice
già esistente, senza scriverne di nuovo:

```matlab
th=sat.geom.parall(i).th;
...
if th>0
    Th=repmat(th,1,6);
    [elem,total_nodes,Connect] = node_box_creator2(Center,L,N,Angles,Th);
else
    [elem,total_nodes,Connect] = node_solid_creator2(Center,L,N,Angles);
end
```

`node_box_creator2.m` è la **stessa identica funzione** già usata per meshare la
Structure esterna (`sat.geom.ext`, con il proprio spessore di parete `th`) — quindi
concettualmente corretto come scelta (un parallelepipedo cavo è geometricamente lo
stesso tipo di oggetto della struttura esterna, solo interno anziché esterno). Il
problema, descritto sotto, è che questo riuso non ha portato con sé una correzione
che esisteva già altrove nel codice proprio per compensare una particolarità di
quella stessa funzione.

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

## 5. Impatto sui risultati esistenti

Qualunque analisi che abbia usato un parallelepipedo **cavo** (`th>0`) avrà uno
scambio radiativo di quell'oggetto dimezzato rispetto a prima — quello vecchio era
sovrastimato del 100%. Trattandosi di una funzionalità introdotta molto di recente
(commit `cb1de4c`, si veda §1), è probabile che non esistano ancora analisi
"validate" storiche che dipendano dal valore precedente.
