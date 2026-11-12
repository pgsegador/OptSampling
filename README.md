# OptSampling

Herramientas de R para optimizar diseños de muestreo de tamaño fijo con
probabilidades de inclusión prescritas. Incluye VNS sobre órdenes sistemáticos
y COS/ACOS mediante símplex revisado y generación de columnas.

Los ejemplos siguientes son sintéticos e independientes de los experimentos del artículo.

## Instalación

Desde la carpeta que contiene el paquete:

```r
install.packages(c("Rcpp", "Matrix"))
install.packages("OptSampling", repos = NULL, type = "source")
library(OptSampling)
```

La instalación compila C++ una sola vez. En Windows hace falta Rtools compatible
con R; en macOS, las herramientas de desarrollo de Xcode; en Linux, un compilador
C++ y las herramientas de desarrollo de R. Requiere R >= 4.1.
También puedes instalarlo directamente desde GitHub:

```r
install.packages("remotes")
remotes::install_github("pgsegador/OptSampling")
```

## Balance con VNS

```r
library(OptSampling)
set.seed(42)
N <- 24L
pik <- rep(6 / N, N)
X <- cbind(size = runif(N, 1, 4), auxiliary = runif(N, 2, 8))
fit <- vns_systematic(pik, X, objective = "balance",
                      n_starts = 3, seed = 17, verbose = FALSE)
print(fit)
summary(fit)
plot(fit, component = fit$selected_runs)
```

`pik` es un vector de longitud N, con valores en (0,1] y suma entera.
`X` tiene una fila por unidad. Cada orden conserva `pik` y el tamaño muestral.
El VNS busca minimizar la **esperanza** de la pérdida sobre todo el soporte,
no la pérdida de una muestra individual. Es una búsqueda local: no garantiza
el óptimo global.

Por defecto, balance minimiza la suma de los errores cuadráticos relativos de
los estimadores HT de los totales. Los divisores son los valores absolutos de
los totales poblacionales. Si algún total es cero, suministra escalas positivas:

```r
scaled <- vns_systematic(pik, X, objective_control = list(scale = c(10, 20)),
                          seed = 17, verbose = FALSE)
projection <- vns_systematic(pik, X,
  objective_control = list(metric = "mahalanobis"), seed = 17, verbose = FALSE)
```

La segunda opción utiliza el criterio de proyección/Mahalanobis y admite
variables redundantes mediante pseudoinversa. `mahalanobis_tol` controla el
umbral numérico de rango. Para minimizar el MSE del total de una variable
auxiliar conocida `z`, usa balance con una columna y escala 1:

```r
z <- X[, 1]
mse_fit <- vns_systematic(pik, matrix(z, ncol = 1),
  objective_control = list(scale = 1), seed = 17, verbose = FALSE)
```

Así `final_objective` está en unidades de MSE, sin normalización. Si `z` es
una predicción, el criterio corresponde al proxy; el MSE real puede ser distinto.

## Arranques y presupuestos de búsqueda

`vns_systematic(n_starts = ...)` ya admite varios órdenes aleatorios.
Con `start = orden` se usa únicamente ese orden y se ignora `n_starts`.
La función adicional admite una lista de órdenes o una matriz con un orden por fila:

```r
set.seed(3)
starts <- replicate(3, sample.int(N), simplify = FALSE)
multi <- vns_multistart(pik, X, starts = starts,
                        time_budgets = c(1, 3), seed = 8, verbose = FALSE)
multi$multistart$best_start
multi$multistart$attempts
multi$multistart$fits
```

También puedes omitir `starts` y usar `n_starts = 10`.
Cada presupuesto se aplica a una llamada de continuación desde el último orden;
**no es un límite acumulado**. Se detiene un arranque cuando obtiene certificado
tipo II o agota la secuencia. Se elige el menor objetivo final; los empates
conservan el primer arranque. La ejecución es secuencial.
`...` permite pasar, por ejemplo, `objective`, `objective_control`, `search`,
`max_iter` e `improvement_tol` a VNS. La estratificación se configura directamente
con `vns_systematic(stratification = TRUE, n_strata = ...)` o `strata = ...`.

`search = "type1"` explora intercambios adyacentes circulares;
`"type2"` explora todas las transposiciones; `"adaptive"` combina ambas búsquedas.
Los certificados son locales y dependen de `improvement_tol`.
`max_iter` limita movimientos aceptados y `max_time` el tiempo por ejecución.
El control de tiempo se realiza periódicamente y puede sobrepasar el presupuesto.
Con límites de tiempo, una semilla no garantiza el mismo punto de parada entre máquinas.

## Leer el resultado y seleccionar una muestra

- `sigma`: orden final de etiquetas 1,...,N.
- `initial_objective`, `final_objective`, `improvement`: criterio e incremento de calidad.
- `diagnostics`: una fila por arranque de la llamada VNS; incluye certificados,
  motivo de parada, tiempo y evaluaciones.
- `runs`: resultados del motor, incluido `f_path`; `selected_runs`: índice ganador.
- En `vns_multistart()`, `multistart$attempts` recoge todas las continuaciones
  y `multistart$fits` el resultado final de cada arranque. Los campos superiores
  `diagnostics` y `runs` describen la última llamada del arranque ganador.

```r
support <- systematic_support(fit$pik, fit$sigma)
support$samples        # lista de vectores de índices de unidades
support$probabilities # probabilidad de cada elemento del soporte
s <- draw_systematic_sample(fit, seed = 99)
X[s, , drop = FALSE]
```

La optimización devuelve un diseño; la extracción aleatoria produce la muestra
que se observaría. Las unidades con `pik = 1` se seleccionan siempre.

## Dispersión espacial

```r
coords <- cbind(runif(N), runif(N))
spatial <- vns_systematic(pik, coords, objective = "spread",
  objective_control = list(metric = "sb"), seed = 4, verbose = FALSE)
s <- draw_systematic_sample(spatial, seed = 5)
plot(coords)
points(coords[s, , drop = FALSE], col = "red", pch = 19)
```

`sb` usa balance espacial de Voronoi; `ib` usa el criterio espacial basado en
Moran. Para `ib`, proporciona `objective_control = list(metric = "ib", W = W)`
con una matriz de proximidad N x N (densa o dispersa). Sin `W`, se requiere
el paquete opcional WaveSampling y se utiliza `WaveSampling::wpik()`.
Las coordenadas deben representar distancias euclídeas adecuadas; para datos
geográficos, proyecta previamente las coordenadas si corresponde.
`dense_distance_limit` (2000 por defecto) y `state_cell_limit` (5e6) permiten
controlar las rutas de cálculo y su consumo de memoria.
`tsp_order(D)` ofrece un orden inicial desde una matriz de distancias y requiere TSP.

## COS y ACOS

```r
pi_cos <- c(.18, .31, .43, .52, .67, .89) # suma 3; todas interiores
values <- c(2, 5, 3, 8, 7, 11)
cos <- cos_design(values, pi_cos, pricing = "exact", seed = 12, verbose = FALSE)
summary(cos)
evaluate_cos_design(cos, values)
draw_cos_sample(cos, seed = 13)
proxy <- values + c(.2, -.3, .4, -.1, .2, -.4)
acos <- acos_design(proxy, pi_cos, seed = 12, verbose = FALSE)
evaluate_cos_design(acos, values)
```

COS minimiza el MSE del estimador HT del total. ACOS llama al mismo motor con
valores predichos en la fase de diseño. `evaluate_cos_design()` evalúa exactamente
el MSE sobre el soporte devuelto. `support` contiene una muestra por fila y
`probabilities` sus probabilidades. `trace` recoge los pivotes; `cos_diagnostics()`
resume factibilidad, objetivo y parada.

`pricing = "exact"` enumera muestras: úsalo sólo en poblaciones pequeñas.
`"heuristic"` usa intercambios; `"hybrid"` decide según `exact_threshold`.
Los controles principales son `max_iter`, `pricing_starts`, `pricing_max_swaps`,
`refactor_every` y las tolerancias (consulta `?cos_design`). Un objetivo evaluado
exactamente sobre el soporte no implica un óptimo global: revisa `certified` y
`status`, especialmente con pricing heurístico o límite de iteraciones.

COS exige 0 < pik < 1. Retira previamente las unidades de certeza, ajusta el
tamaño residual y añádelas después a cada muestra; su contribución al total es
determinista. VNS admite directamente probabilidades de certeza.

## Otras herramientas

`estimate_sampling_design_mc()` estima momentos de un generador de muestras;
`estimate_cube_design()` lo hace para cube y requiere sampling.
`approximate_design_by_systematics()` aproxima los momentos de segundo orden
mediante una mezcla de diseños sistemáticos. Al pasar esa mezcla a VNS, sus
componentes se usan como órdenes iniciales y el resultado final es un solo diseño.
`cos_pps_inclusion_probabilities()` construye probabilidades PPS y
`cos_correlated_proxy()` genera un proxy con correlación prescrita.

## Estado de publicación

Versión inicial de desarrollo. La licencia de redistribución está pendiente
de elección por los titulares del código; consulta LICENSE.
No se incluyen el artículo, sus datos, resultados ni scripts de experimentos.
