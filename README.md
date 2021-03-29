# roots
Subrutinas en Fortran para la resolución de ecuaciones no lineales de una variable

En general, las raíces de una ecuación no lineal f (x) = 0 no pueden ser obtenidas por fórmulas explícitas cerradas, con lo que no es posible obtenerlas en forma exacta. De este modo, para resolver la ecuación nos vemos obligados a obtener soluciones aproximadas a través de algún método numérico.

En el módulo Fortran roots implementamos los métodos numéricos usuales, como ser bisección, Newton, secante y punto fijo, junto a un método de propósito general más robusto y eficiente.
