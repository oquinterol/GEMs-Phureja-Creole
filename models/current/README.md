# Modelo Actual: creole_v1.9_consistency_fixed.xml

Este directorio contiene la versión más reciente y en uso del modelo metabólico, utilizada para las simulaciones actuales.

## Versión: 1.9 (Consistencia Corregida)

Esta versión representa una actualización crítica con respecto a la `v1.8`, enfocada principalmente en la corrección de la consistencia interna del modelo para mejorar la fiabilidad de las simulaciones.

### Resumen de Cambios:

El nombre del archivo, `consistency_fixed`, refleja los cambios principales realizados:

1.  **Corrección de la Direccionalidad de Reacciones**: Se ha llevado a cabo una revisión de las reacciones de intercambio (`Exchange Reactions`). Se corrigió la direccionalidad de varias de ellas, cambiando la definición de `listOfReactants` a `listOfProducts` para asegurar que el modelo pueda importar metabolitos desde el medio externo de forma correcta durante las simulaciones.

2.  **Simplificación de Anotaciones**: Se han simplificado las etiquetas de los productos génicos (`geneProduct`). Las anotaciones detalladas (KEGG, GO terms) han sido eliminadas de las etiquetas, lo que resulta en un archivo de modelo más limpio, ligero y fácil de procesar por las herramientas de simulación.

3.  **Curación del Modelo**: Se han eliminado algunas especies y reacciones específicas que estaban presentes en la v1.8 como parte de un proceso de curación para mejorar la calidad y consistencia general del modelo.

En resumen, la **versión 1.9** es una versión más robusta y consistente, optimizada para análisis computacionales como FBA (Flux Balance Analysis).