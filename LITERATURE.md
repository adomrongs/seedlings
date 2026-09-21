# Literatura

Ordenada por cercanía a la pregunta del proyecto.

## 1. El núcleo: optimal cross / mate selection y el sesgo contra los no fenotipados

**Endelman (2025), *Genomic prediction of heterosis, inbreeding control, and mate
allocation in outbred diploid and tetraploid populations*, Genetics 229(2):iyae193.**
El paper del paquete COMA que ya usas. Es el punto de partida obligado: define la
optimización convexa sobre cruces, la restricción de dF y la parametrización
'Breeding'. Lo relevante aquí es que el mérito del cruce se construye a partir de
los méritos parentales predichos, que es exactamente donde entra el shrinkage
diferencial que quieres medir.
<https://academic.oup.com/genetics/article/229/2/iyae193/7903017>

**Gorjanc, Gaynor & Hickey (2018), *Optimal cross selection for long-term genetic
gain in two-part programs with rapid recurrent genomic selection*, TAG 131:1953-1966.**
El precedente más directo de tu pregunta. Muestran que la optimización de cruces
sostiene la ganancia a largo plazo cuando los candidatos son material joven
seleccionado genómicamente. Su "two-part" es justamente un programa donde el
population improvement usa individuos sin fenotipo.
<https://link.springer.com/article/10.1007/s00122-018-3125-3>

**Allier et al. (2019), *Improving short and long term genetic gain by accounting
for within-family variance in optimal cross selection*, Front. Genet. 10:1006.**
Añaden la varianza intra-familia al criterio. Importante para ti porque la
varianza predicha de una familia de seedlings no sufre el mismo shrinkage que la
media, así que puede ser parte de la respuesta a tu pregunta 2.
<https://www.frontiersin.org/articles/10.3389/fgene.2019.01006/full>

**Woolliams, Berg, Dagnachew & Meuwissen (2015), *Genetic contributions and their
optimization*, J. Anim. Breed. Genet. 132:89-99**, y **Clark et al. (2013), *The
effect of genomic information on optimal contribution selection in livestock
breeding programs*, GSE 45:44.**
El segundo es el que te interesa más: muestra cómo cambia OCS cuando la precisión
de la predicción es heterogénea entre candidatos. Es el análogo en ganado de tu
pregunta.
<https://link.springer.com/article/10.1186/1297-9686-45-44>

**Christensen et al. (2025), *Uncertainty-aware breeding decisions: MCMC-based
optimum contribution selection increases breeding decision robustness*, Genetics.**
Muy pertinente: trata explícitamente la incertidumbre de los EBVs dentro de OCS
en lugar de tratarlos como valores fijos. Es una vía alternativa a la de fenotipar
seedlings para corregir el sesgo.
<https://doi.org/10.1093/genetics/iyag205>

## 2. GS en programas clonales, que es tu esquema

**Wu, Chen, Stich et al. (2023, 2024), *Optimal implementation of genomic selection
in clone breeding programs — exemplified in potato*, The Plant Genome, partes I y II.**
La referencia más cercana a tu diseño. Parte I: en qué etapa implementar GS y con
qué intensidad de selección, efecto sobre la ganancia a corto plazo. Parte II:
estrategia de selección y método de selección de cruces sobre la ganancia a largo
plazo. Tu decisión de aplicar GS solo de S2 a S3 debería contrastarse con lo que
encuentran.
Parte I: <https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20327>
Parte II: <https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.70000>

**Gaynor et al. (2017), *A two-part strategy for using genomic selection to develop
inbred lines*, Crop Science 57:2372-2386.**
El paper que popularizó la separación entre population improvement y product
development. Tu pregunta se puede releer como: ¿hasta qué punto un programa clonal
consigue, de facto, esa separación si la optimización de cruces excluye a los
seedlings?
<https://acsess.onlinelibrary.wiley.com/doi/10.2135/cropsci2016.09.0742>

**Grattapaglia et al. (2018), *Expected benefit of genomic selection over forward
selection in conifer breeding and deployment*, PLOS ONE 13(12):e0208232.**
Cuantifica la ganancia por año al acortar el ciclo. Útil para poner en contexto
cuánto cuesta que los seedlings no entren en el mating plan.
<https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208232>

## 3. Herramientas y simulación

**Gaynor, Gorjanc & Hickey (2021), *AlphaSimR: an R package for breeding program
simulations*, G3 11(2):jkaa017.** La cita del motor.
<https://academic.oup.com/g3journal/article/11/2/jkaa017/6025179>

**Bančič et al. (2025), *Plant breeding simulations with AlphaSimR*, Crop Science.**
Tutorial reciente y muy práctico sobre cómo estructurar burnin, escenarios futuros
y comparaciones. Buen espejo para validar tus decisiones de diseño.
<https://acsess.onlinelibrary.wiley.com/doi/10.1002/csc2.21312>

**Peixoto et al. (2025), *SimpleMating: R-package for prediction and optimization of
breeding crosses using genomic selection*, The Plant Genome.**
Alternativa a COMA con un criterio de usefulness distinto. Vale la pena mirarlo
aunque no lo uses, para justificar por qué COMA.
<https://acsess.onlinelibrary.wiley.com/doi/10.1002/tpg2.20533>

**COMA, repositorio.** <https://github.com/jendelman/COMA>

## 4. Diseño de training population y precisión

**Lorenz & Smith (2015) y Rio et al. (2019)** sobre el efecto del parentesco entre
training y candidatos en la precisión. Es el mecanismo detrás del shrinkage
diferencial que quieres medir: tus S2 son, por construcción, la generación más
alejada de una TP hecha de S3 y S4 de hasta 8 años atrás.
Ver también Werner et al. (2020) sobre cómo la estructura de la población afecta
a la precisión estimada por validación cruzada.
<https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2020.592977/full>

---

## Lectura sugerida para empezar

Endelman 2025 (mecanismo) → Wu 2023/2024 (tu esquema) → Gorjanc 2018 (tu pregunta
en otro contexto) → Christensen 2025 (la alternativa que no habías considerado).
