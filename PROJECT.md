# seedlings

Simulación de un programa de mejora clonal para estudiar **qué proporción de
seedlings no fenotipados entra en el mating plan** cuando la candidate pool
mezcla individuos genotipados de varias etapas del ciclo.

Autor: Alejandro Domínguez · TheRocinanteLab · septiembre 2026
Proyecto hermano: `OMA_Comparison` (autotetraploide, COMA con dominancia).

---

## 1. Pregunta

Cuando la candidate pool mezcla individuos con y sin fenotipo, el modelo aplica
más *shrinkage* a los no evaluados. Sus GEBVs quedan empujados hacia la media
poblacional y, en consecuencia, también las medias de las familias generadas a
partir de ellos. Esto sesga el optimum mating plan contra los seedlings,
justo los individuos que harían corto el ciclo de mejora. Las preguntas del
documento `Crossing_Designs.pdf`:

1. ¿En qué rango se mueve la proporción de individuos no fenotipados incluidos
   en un optimum mating plan?
2. Si esa proporción es siempre baja, ¿un mating plan con más seedlings es
   mejor que el optimum mating plan?
3. Si no lo es, ¿fenotipar una fracción de los seedlings aumenta su presencia
   en el mating plan, y aumenta eso la ganancia?

**Este repositorio aborda la pregunta 1.** Las preguntas 2 y 3 requieren
escenarios adicionales que todavía no están implementados (ver §8).

---

## 2. Esquema del programa

```
Crossing block   40 padres  ->  100 cruces x 200 seedlings
S1               20 000 seedlings/año      h2 = 0.1
S2                1 000 individuos         h2 = 0.2
S3                  200 clones             h2 = 0.4
S4                   20 clones             h2 = 0.8
```

Una etapa por año. El pipeline queda desfasado: F1 es la cohorte más reciente y
S4 la más antigua, con cuatro años entre el cruce y la evaluación final.

- **Burnin**: 20 años, **todas las selecciones fenotípicas**. El crossing block
  se recicla sustituyendo los 10 padres más antiguos por los 10 mejores S4.
- **Futuro**: 30 años, 30 repeticiones. Cada repetición parte del mismo burnin.

### Heredabilidades

Se fijan con `setPheno(h2 = ...)`. AlphaSimR deriva `varE` de la varianza
genética del founder pop, así que `varE` es constante y la h2 realizada baja
con los años a medida que se agota la varianza genética. Es el comportamiento
deseado: lo que se mantiene fijo es el error experimental, no la heredabilidad.

Un efecto de año común a todas las etapas se introduce con el argumento `p` de
`setPheno()`, muestreado una sola vez por repetición (`P <- runif(50)`).

---

## 3. Escenarios

| Escenario | Selección S2→S3 | Crossing block | Script |
|---|---|---|---|
| `burnin` | fenotípica | 10 mejores S4 por fenotipo | `code/scenarios/burnin.R` |
| `PS` | fenotípica | 10 mejores S4 por fenotipo | `code/scenarios/PS.R` |
| `GS_trunc` | **GEBV** | 40 mejores GEBV de S2+S3+S4 | `code/scenarios/GS_trunc.R` |
| `GS_COMA_0.5` | **GEBV** | COMA `oma()` sobre S2+S3+S4, dF = 0.5 % | `code/scenarios/GS_COMA_0.5.R` |

`PS` no estaba en el encargo original: es la continuación del burnin durante los
mismos 30 años y sirve de línea base para medir la ganancia de los dos
escenarios GS. Si no se quiere, basta con borrar el script y su `.slurm`.

En los tres escenarios futuros **solo S2→S3 es genómica**. S1→S2 y S3→S4 siguen
siendo fenotípicas, tal como se especificó.

---

## 4. Decisiones cerradas

| Decisión | Valor | Motivo |
|---|---|---|
| Ploidía y modelo | Diploide, **solo aditivo** | Código simple. Sin dominancia, sin `Compute_Q`, sin StageWise |
| Genoma | 10 crom., 1.43 M, 8e8 pb, µ = 2e-9, 1500 seg. sites/crom. | Covarrubias et al. 2026, sin líneas puras |
| Marcadores | 100 QTL + 1000 SNP por cromosoma | 1000 QTL y 10 000 SNP en total |
| Training population | **S3 + S4**, ventana de 8 años | 220/año × 8 = 1760 registros, por encima del mínimo de 1500 |
| Candidate population | **Padres del plan anterior + S2 + S3 + S4**, sin filtrar | 1220 de las etapas más los padres que no hayan avanzado. Ver §6 |
| Modelo de predicción | GBLUP: `AGHmatrix::Gmatrix` (VanRaden) + `lme4breeding::lmebreed` | `pheno ~ year + stage + (1 \| gid)` |
| Efectos de marcador | Back-solve de los GEBVs | Necesarios para predecir candidatos fuera del modelo y para alimentar COMA |
| Genotipado | S2, S3 y S4 cada año, solo en escenarios futuros | El burnin no genotipa. Los padres reciclados ya se genotiparon en su día; el crossing block heredado del burnin se genotipa una vez, en el año 21 |
| Inbreeding | Pedigrí (`polyBreedR::A_mat`) | El pedigrí se guarda completo desde el año 1 |
| Ejecución | SLURM, array 1-30, `SLURM_ARRAY_TASK_ID` con fallback a 1 | Corre también en local sin tocar nada |

### Modelo de predicción, en detalle

```
pheno ~ year + stage + (1 | gid),  relmat = list(gid = G)
```

`stage` entra como efecto fijo porque S3 y S4 tienen heredabilidades distintas
(0.4 y 0.8) y, por tanto, medias y precisiones distintas. `year` y `stage` son
estimables por separado porque cada año aporta registros de ambas etapas. Un
mismo clon aparece en S3 un año y en S4 al siguiente, así que el modelo es de
medidas repetidas sobre `gid`.

Los efectos de marcador se recuperan de los GEBVs del training con

```
a = W' (W W' + I·1e-6)^-1 · GEBV
```

donde `W` es la matriz de dosis centrada. Es lo que permite predecir S2 (que no
tiene fenotipo y no está en el modelo) y lo que COMA espera en `geno.file`.

---

## 5. Métricas registradas

Una fila por año y repetición, en `outputs/<escenario>/csv/<escenario>_<rep>.csv`.

| Columna | Definición |
|---|---|
| `meanG_S1`, `varG_S1`, `genicVarA_S1` | Ganancia, varianza genética y **varianza genética aditiva** en S1 |
| `accuracy_S2` | Precisión de selección en la decisión S2→S3: `cor(gv, ebv)` (o `cor(gv, pheno)` en PS) |
| `accuracy_cand` | `cor(gv, ebv)` en toda la candidate pool |
| `Ft0`, `Ft1`, `dF1` | Inbreeding del crossing block, de la F1, y tasa `100·(Ft1−Ft0)/(1−Ft0)` |
| `propParents_S2` | % de padres del mating plan que son S2 |
| `propCrosses_S2` | % de cruces con al menos un padre S2 |
| `propProgeny_S2` | % de la progenie total asignada a padres S2 |
| `nCandidates`, `nParentsUsed`, `nCrossesUsed`, `nKept` | Tamaños y solapamiento del crossing block entre años |
| `time` | Minutos de `COMA::oma()` |

Las tres proporciones responden la pregunta 1 desde ángulos distintos.
`propProgeny_S2` es la más informativa bajo COMA, porque pondera por la
contribución óptima de cada cruce; bajo truncación las tres coinciden salvo por
el redondeo, ya que todos los cruces reciben la misma progenie.

En `burnin` y `PS` estas tres columnas quedan `NA`: por diseño ningún S2 puede
entrar en el crossing block, así que no hay nada que medir.

`Ft1` se calcula como la **kinship media de los pares de padres** del mating
plan, no construyendo la matriz A de los 20 000 seedlings. Es exacto
(`F_hijo = kinship(p1, p2)`) y evita una matriz de 3 GB.

---

## 6. Aviso de coste computacional

La candidate pool son los 1220 individuos de las etapas (1000 + 200 + 20) más
los padres del plan del año anterior que no hayan avanzado de etapa. Ese
segundo término no es despreciable y **no es constante**: en el año 21 entra el
crossing block del burnin al completo (40), y en años posteriores entran los
padres que no pasaron a S3 o S4. Con COMA, que puede usar bastantes más de 40
padres, el añadido es mayor y varía año a año.

Los cruces escalan con el cuadrado de la pool, así que conviene tenerlo medido:

| Pool | Cruces tras quitar recíprocos y selfs | Factor |
|---|---|---|
| 1220 (solo etapas) | 743 590 | 1.00× |
| 1260 (+40, típico de `GS_trunc`) | 793 170 | 1.07× |
| 1320 (+100) | 870 540 | 1.17× |
| 1420 (+200, COMA con plan amplio) | 1 007 490 | 1.35× |

`nS2` se fijó en 1000, y no en 2000, precisamente por esto: con 2000 la pool de
etapas sube a 2220 y los cruces a 2 463 090, **3.3 veces más**, y COMA deja de
ser ejecutable a esta escala.

La columna `nCandidates` del CSV registra el tamaño real de la pool cada año.
Si crece de forma sostenida porque COMA retiene muchos padres, el coste sube
con el cuadrado: es la primera señal que hay que mirar si los trabajos empiezan
a no caber en el tiempo pedido.

**Antes de lanzar el array completo, ejecutar una sola repetición un solo año y
mirar la columna `time`.** Si un año tarda más de ~30 minutos, hay que decidir
entre: prefiltrar S2 por GEBV, bajar `nS2` otra vez, o aceptar el coste. El
`--time=48:00:00` y los `--mem=120G` del `.slurm` son estimaciones sin dato
empírico detrás, dimensionadas para el caso de 2220 candidatos: pedir de más
solo cuesta tiempo de cola, así que se dejan hasta tener el benchmark.

---

## 7. Estructura y ejecución

```
seedlings/
├── seedlings.Rproj
├── PROJECT.md          <- este documento
├── PROMPT.md           <- prompt para retomar o extender el proyecto
├── LITERATURE.md
├── code/
│   ├── setup_renv.R
│   ├── processes/      <- una acción por script, se sourcean desde los escenarios
│   ├── scenarios/      <- un script por escenario, es lo que se ejecuta
│   └── slurm/
└── outputs/            <- ignorado por git salvo .gitkeep
```

### Procesos

| Script | Acción |
|---|---|
| `0_params.R` | Parámetros globales y cuatro funciones auxiliares |
| `1_CreateFounders.R` | Founder pop, rasgo aditivo, chip de SNPs, crossing block inicial |
| `2_FillPipeline.R` | Rellena el pipeline con 5 cohortes desfasadas |
| `3_UpdateParents.R` | Recicla el crossing block por fenotipo |
| `4_AdvanceYear.R` | Avanza un año, todo fenotípico |
| `5_StoreRecords.R` | Registros de S3+S4 (ventana de 8 años) y pedigrí acumulado |
| `6_Inbreeding.R` | `Ft0`, `Ft1`, `dF1` |
| `7_PredictGEBV.R` | GBLUP, efectos de marcador, GEBVs de la candidate pool (padres + S2 + S3 + S4) |
| `8_SelectParentsTrunc.R` | Crossing block por truncación sobre GEBV |
| `9_COMAFiles.R` | `geno.file` y `kinship.file` para COMA |
| `10_RunOMA.R` | `COMA::oma()` y `SimPlus::sim_mate()` |
| `11_AdvanceYearGS.R` | Avanza un año con GS de S2 a S3 |

Los procesos operan sobre objetos del entorno global, como en `OMA_Comparison`.
No son funciones: son bloques de acción que los escenarios encadenan con
`source()`. Las únicas funciones están en `0_params.R`.

### Orden dentro de un año

```
burnin / PS         3_UpdateParents -> 4_AdvanceYear -> 5_StoreRecords -> 6_Inbreeding
GS_trunc            7_PredictGEBV -> 8_SelectParentsTrunc -> 11_AdvanceYearGS -> 5_StoreRecords -> 6_Inbreeding
GS_COMA_0.5         7_PredictGEBV -> 9_COMAFiles -> 10_RunOMA -> 11_AdvanceYearGS -> 5_StoreRecords -> 6_Inbreeding
```

`6_Inbreeding` va siempre después de `5_StoreRecords`, porque es este el que
mete a los padres del año en el fichero de pedigrí.

### Ejecutar

Siempre desde la raíz del proyecto (`seedlings.Rproj` abierto o `setwd()`).

```bash
sbatch code/slurm/burnin.slurm          # primero, genera los .Rdata de partida
sbatch code/slurm/PS.slurm
sbatch code/slurm/GS_trunc.slurm
sbatch code/slurm/GS_COMA_0.5.slurm
```

En local, sin SLURM, `rep` cae a 1:

```r
source('code/scenarios/burnin.R')
```

### renv

`code/setup_renv.R` se ejecuta **una sola vez**, a mano. Usa `renv::init(bare = TRUE)`
en lugar de `renv::init()` a secas: el init normal escanea el código, encuentra
los `library()` de paquetes que solo están en GitHub, falla al buscarlos en CRAN
y deja el proyecto a medias. Instalando explícitamente con `renv::install('user/repo')`
el remote queda anotado en `renv.lock` y `renv::restore()` funciona en el cluster.

**Pendiente**: el remote de `SimPlus` está como TODO en `setup_renv.R`. Solo lo
usa `10_RunOMA.R`.

---

## 8. Qué falta

- Escenario con mating plan forzado a más seedlings (pregunta 2 del PDF).
- Escenario donde se fenotipa una fracción de S1/S2 y entra en la training
  population (pregunta 3).
- Escenario de sensibilidad al número total de seedlings disponibles.
- Script de agregación y figuras (`code/results.R`).

---

## 9. Coding guidelines

Guías de comportamiento para reducir errores habituales. Se combinan con las
instrucciones específicas del proyecto.

*Compromiso: estas guías priorizan prudencia sobre velocidad. Para tareas
triviales, usa el criterio.*

### 9.1 Pensar antes de programar

No asumas. No escondas la confusión. Saca los compromisos a la luz.

Antes de implementar:

- Enuncia tus supuestos de forma explícita. Si hay incertidumbre, pregunta.
- Si existen varias interpretaciones, preséntalas: no elijas en silencio.
- Si existe un enfoque más simple, dilo. Discrepa cuando esté justificado.
- Si algo no está claro, para. Nombra qué te confunde. Pregunta.

### 9.2 Simplicidad primero

El mínimo código que resuelve el problema. Nada especulativo.

- Ninguna funcionalidad más allá de lo pedido.
- Ninguna abstracción para código de un solo uso.
- Ninguna "flexibilidad" o "configurabilidad" que no se haya pedido.
- Ningún manejo de errores para escenarios imposibles.
- Si escribes 200 líneas y podrían ser 50, reescríbelo.

Pregúntate: "¿un ingeniero senior diría que esto está sobrecomplicado?" Si sí,
simplifica.

### 9.3 Cambios quirúrgicos

Toca solo lo imprescindible. Limpia solo tu propio desorden.

Al editar código existente:

- No "mejores" código, comentarios o formato adyacentes.
- No refactorices lo que no está roto.
- Respeta el estilo existente, aunque tú lo harías de otra forma.
- Si detectas código muerto no relacionado, menciónalo: no lo borres.

Cuando tus cambios dejan huérfanos:

- Elimina imports, variables y funciones que **tus** cambios dejaron sin uso.
- No elimines código muerto preexistente salvo que te lo pidan.

La prueba: cada línea cambiada debe poder trazarse directamente a lo que se pidió.

### 9.4 Ejecución dirigida por objetivos

Define criterios de éxito. Itera hasta verificarlos.

Convierte tareas en objetivos verificables:

- "Añade validación" → "Escribe tests para entradas inválidas y haz que pasen"
- "Arregla el bug" → "Escribe un test que lo reproduzca y haz que pase"
- "Refactoriza X" → "Asegura que los tests pasan antes y después"

Para tareas de varios pasos, enuncia un plan breve:

```
1. [Paso] → verificar: [comprobación]
2. [Paso] → verificar: [comprobación]
3. [Paso] → verificar: [comprobación]
```

Criterios de éxito fuertes permiten iterar de forma autónoma. Criterios débiles
("haz que funcione") obligan a aclarar constantemente.

### 9.5 Verificar que el cambio ha llegado

Nunca des un cambio por hecho sin comprobarlo. "He editado el fichero" no es lo
mismo que "el fichero de destino contiene el cambio".

Antes de decir que algo está implementado:

- **Lee de vuelta el resultado desde su destino final**, no desde la copia
  intermedia que acabas de escribir. Si el fichero viaja entre máquinas, la
  única fuente de verdad es el fichero en destino.
- **Compara**: checksum contra el original, o `grep` de las líneas concretas que
  debían cambiar. Un vistazo al fichero que tú mismo escribiste no verifica nada.
- **Un código de éxito no es una verificación.** Una herramienta puede devolver
  "escrito, sin errores" y haber transferido bytes obsoletos. Ya ha pasado en
  este proyecto: una transferencia reportó éxito y dejó la versión anterior en
  disco. Trata la respuesta de la herramienta como un acuse de recibo, no como
  una confirmación del contenido.
- **Riesgo alto de caché al reescribir.** Si el fichero ya se había transferido
  antes a esa misma ruta, la probabilidad de recibir una copia obsoleta sube.
  Ante la duda, reintenta desde una ruta intermedia nueva y vuelve a verificar.
- **Si la verificación falla, dilo y arréglalo**, no lo silencies ni lo reportes
  como hecho.

Lo mismo aplica al código: un script que pasa `parse()` no está verificado. Hay
que comprobar que los objetos que consume existen en el orden en que se sourcea
dentro del escenario, con el arnés de prueba o con una ejecución reducida.

**Estas guías funcionan si**: hay menos cambios innecesarios en los diffs, menos
reescrituras por sobrecomplicación, las preguntas de aclaración llegan antes de
implementar en lugar de después de equivocarse, y ningún cambio se reporta como
hecho sin haberlo leído en su destino.
