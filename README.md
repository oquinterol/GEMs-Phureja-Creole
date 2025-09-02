# Reconstrucción y Validación de GEMs — *Solanum tuberosum* Gr. *Phureja* (var. Colombia)

Repositorio para la **reconstrucción, curación y validación** de Modelos Metabólicos a Escala Genómica (GEMs) de la variedad **Criolla Colombia**. El proyecto integra **COBRApy**, **ModelSEEDpy**, y **MEMOTE**, con cumplimiento estricto de **SBML Level 3 V2 + FBC v2** y **anotaciones MIRIAM/Identifiers.org**. Está diseñado para ejecutarse **localmente** (offline) o mediante **APIs** (resolución y enriquecimiento de metadatos en línea).

---

## Objetivos

1. **Reconstruir** un GEM de *S. tuberosum* Gr. *Phureja* (Cultivar Criolla Colombia) con soporte genómico.
2. **Curar** y **anotar** el modelo con MIRIAM (Identifiers.org), términos SBO y referencias externas (ChEBI/Rhea/KEGG/UniProt/GO).
3. **Validar** calidad y reproducibilidad con **MEMOTE** y *smoke tests* en **COBRApy**.
4. **Publicar** artefactos y *reports* trazables (FAIR), con opción de archivo **COMBINE (OMEX)** y **CITATION.cff**.

---

## Características

- **SBML L3 V2 + FBC v2**: restricciones de flujo, objetivos, GPRs y *bounds* explícitos.
- **MIRIAM/Identifiers.org**: URIs persistentes para metabolitos (ChEBI), reacciones (Rhea/KEGG), genes/proteínas (UniProt) y compartimentos (GO).
- **MEMOTE**: informe HTML/JSON y umbral recomendado ≥ 80.
- **COBRApy**: pruebas de humo de FBA (solver GLPK por defecto).
- **Ejecución dual**: *Local* (offline) o *API-mode* (resolución de metadatos).


## Convenciones y estándares

- **SBML L3V2 + FBC v2** para el modelo final.
- **MIRIAM/Identifiers.org** para URIs persistentes (ChEBI, Rhea/KEGG, UniProt, GO).
- **SBO** añadido a reacciones, especies y objetivo.
- **Nomenclatura**:
  - Metabolitos: `id_comp` (p. ej., `glc__D_c`).
  - Compartimentos: `{c,m,p,x,v,e,g,r}`.
  - Reacciones: `R_*` con SBO y referencias externas.
  - Genes/Proteínas: `geneProduct` (FBC) con mapeo UniProt/Ensembl si aplica.

---

## Contribuir

1. Abre un *issue* describiendo el cambio propuesto.
2. Crea una rama desde `main` y realiza *commits* atómicos.
3. Incluye pruebas (smoke/MEMOTE) y actualiza documentación si aplica.
4. Envía un *pull request* con un resumen claro del cambio y su justificación.

> Se espera un estilo de código claro, notebooks con propósito definido y resultados reproducibles. Por favor, evita subir datos con restricciones de licencia; comparte únicamente identificadores y metadatos mínimos cuando sea necesario.

---

## Agradecimientos

A la comunidad de **COBRA**, **MEMOTE**, **ModelSEED/ModelSEEDpy**, **CarveMe**, y a los consorcios que mantienen recursos como **UniProt**, **Rhea**, **ChEBI**, **GO**, **KEGG**, y **PMN/PlantCyc**. El uso de ciertas APIs/datos puede estar sujeto a términos específicos; por favor, revisa y respeta las licencias asociadas.
