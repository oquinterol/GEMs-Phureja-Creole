# STRING Protein-Protein Interaction Network Analysis

Este directorio contiene scripts para analizar redes de interacciones proteína-proteína (PPI) usando la base de datos STRING para *Solanum tuberosum* Group Phureja.

## Contexto del Problema

**Problema Original**: Intentaste cargar `Predicted_Proteins.faa` (~205K secuencias) directamente a la interfaz web de STRING, pero recibiste un error porque:

1. **Límite de STRING**: Máximo 2,000 proteínas por upload en la interfaz web
2. **Tu archivo**: 204,709 secuencias de proteínas
3. **Comportamiento de STRING**: Colapsa múltiples isoformas al nivel de gen (una proteína canónica por locus)

## Solución Implementada

Hemos creado un pipeline de dos pasos:

### 1. Filtrado Inteligente del Proteoma (`filter_proteome_for_string.py`)

Este script:
- ✅ Filtra proteínas por relevancia metabólica usando anotaciones de eggNOG-mapper
- ✅ Selecciona la isoforma más larga por gen (mimicking el comportamiento de STRING)
- ✅ Divide las proteínas en lotes de tamaño configurable (máx 2000 para web)
- ✅ Genera reportes de resumen con estadísticas de anotaciones

### 2. Análisis Automatizado via API (`string_api_batch.py`)

Este script:
- ✅ Consulta la API REST de STRING de forma programática
- ✅ Mapea identificadores de proteínas a STRING IDs
- ✅ Recupera redes PPI con scores de confianza
- ✅ Obtiene análisis de enriquecimiento funcional (GO, KEGG, etc.)
- ✅ Exporta resultados en formato TSV

---

## Instalación de Dependencias

```bash
# Activar el ambiente pyenv
source ~/.pyenv/versions/Creole/bin/activate

# Instalar requests para el script de API
pip install requests
```

---

## Uso Rápido

### Opción 1: Exploración Manual (Interfaz Web de STRING)

**Mejor para**: Exploración inicial, visualización interactiva de redes pequeñas

```bash
cd /home/oscar/Documents/Master/Thesis/GEM_Creole

# Los archivos ya están generados en data/string/web_example/
# Usa creole_web_example_batch01.faa para empezar (200 proteínas con EC numbers)

# Pasos:
# 1. Ve a https://string-db.org/
# 2. Click en "Multiple proteins"
# 3. Selecciona organismo: "Solanum tuberosum" (taxonomy ID: 4113)
# 4. Upload: data/string/web_example/creole_web_example_batch01.faa
# 5. Explora la red interactivamente
```

**Archivos disponibles para web**:
- `data/string/web_example/creole_web_example_batch01.faa` - 200 proteínas (RECOMENDADO para empezar)
- Archivos batch02-60 con 200 proteínas cada uno
- Total: 11,832 proteínas con EC numbers (alta confianza)

### Opción 2: Análisis Programático (API de STRING)

**Mejor para**: Análisis sistemático a gran escala, integración con pipelines

```bash
cd /home/oscar/Documents/Master/Thesis/GEM_Creole

# Analizar el primer batch de ejemplo
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch01.faa \
  --out data/string/results \
  --species 4113 \
  --min-score 400 \
  --get-enrichment

# Resultados generados:
# - creole_web_example_batch01_network.tsv      # Interacciones PPI
# - creole_web_example_batch01_enrichment.tsv   # Enriquecimiento funcional
# - creole_web_example_batch01_string_ids.tsv   # Mapeo de IDs
```

**Analizar múltiples batches en secuencia**:

```bash
# Procesar los primeros 5 batches
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch0{1,2,3,4,5}.faa \
  --out data/string/results \
  --species 4113 \
  --min-score 700  # Alta confianza
```

---

## Guía Detallada de Scripts

### 🔧 Script 1: `filter_proteome_for_string.py`

#### Descripción
Prepara el proteoma para análisis en STRING mediante filtrado inteligente y batching.

#### Sintaxis Básica

```bash
python scripts/string/filter_proteome_for_string.py \
  --proteome <FASTA_FILE> \
  --annotations <EMAPPER_FILE> \
  --out <OUTPUT_DIR> \
  [--batch-size SIZE] \
  [--min-score SCORE] \
  [--require-ec] \
  [--require-kegg] \
  [--prefix PREFIX]
```

#### Parámetros

| Parámetro | Tipo | Descripción | Default |
|-----------|------|-------------|---------|
| `--proteome` | Path | Archivo FASTA del proteoma | **Requerido** |
| `--annotations` | Path | Archivo de anotaciones eggNOG-mapper | **Requerido** |
| `--out` | Path | Directorio de salida para batches | **Requerido** |
| `--batch-size` | Int | Máximo de secuencias por archivo | 2000 |
| `--min-score` | Float | Score mínimo de anotación (0-1000) | 50.0 |
| `--require-ec` | Flag | Requiere número EC | False |
| `--require-kegg` | Flag | Requiere anotación KEGG | False |
| `--prefix` | String | Prefijo para archivos de salida | creole_metabolic |

#### Ejemplos de Uso

```bash
# 1. Filtro permisivo: Cualquier anotación KEGG
python scripts/string/filter_proteome_for_string.py \
  --proteome data/Proteome/Predicted_Proteins.faa \
  --annotations data/annotations/emapper_creole.emapper.annotations \
  --out data/string/batches \
  --require-kegg \
  --batch-size 500
# Resultado: ~34,695 proteínas en 70 batches

# 2. Filtro estricto: Solo proteínas con EC numbers y score alto
python scripts/string/filter_proteome_for_string.py \
  --proteome data/Proteome/Predicted_Proteins.faa \
  --annotations data/annotations/emapper_creole.emapper.annotations \
  --out data/string/high_confidence \
  --require-ec \
  --min-score 100 \
  --batch-size 200
# Resultado: ~11,832 proteínas en 60 batches

# 3. Filtro para web STRING (lotes pequeños para exploración)
python scripts/string/filter_proteome_for_string.py \
  --proteome data/Proteome/Predicted_Proteins.faa \
  --annotations data/annotations/emapper_creole.emapper.annotations \
  --out data/string/web_batches \
  --require-ec \
  --batch-size 100 \
  --prefix web_ready
# Resultado: Lotes de 100 proteínas (ideal para exploración web)
```

#### Archivos de Salida

```
<OUTPUT_DIR>/
├── <prefix>_batch01.faa         # Primer lote de secuencias
├── <prefix>_batch02.faa         # Segundo lote
├── ...
├── <prefix>_batchNN.faa         # Último lote
└── <prefix>_summary.txt         # Reporte de resumen
```

**Contenido del reporte de resumen**:
- Total de secuencias originales
- Genes únicos identificados
- Transcripts que pasaron filtros
- Cobertura de anotaciones (EC, KEGG, GO)
- Estrategia de filtrado aplicada

---

### 🌐 Script 2: `string_api_batch.py`

#### Descripción
Consulta la API REST de STRING para obtener redes PPI y análisis de enriquecimiento funcional.

#### Sintaxis Básica

```bash
python scripts/string/string_api_batch.py \
  --input <FASTA_FILE(S)> \
  --out <OUTPUT_DIR> \
  [--species TAXON_ID] \
  [--min-score SCORE] \
  [--get-enrichment] \
  [--network-type TYPE]
```

#### Parámetros

| Parámetro | Tipo | Descripción | Default |
|-----------|------|-------------|---------|
| `--input` | Path(s) | Archivo(s) FASTA a procesar (acepta múltiples) | **Requerido** |
| `--out` | Path | Directorio de salida | **Requerido** |
| `--species` | Int | NCBI Taxonomy ID | 4113 (*S. tuberosum*) |
| `--min-score` | Int | Score mínimo de confianza (0-1000) | 400 (medium) |
| `--get-enrichment` | Flag | Obtener análisis de enriquecimiento | False |
| `--network-type` | String | Tipo de red: functional o physical | functional |

#### Scores de Confianza STRING

| Score | Nivel | Descripción |
|-------|-------|-------------|
| 150-400 | Low | Confianza baja, muchas predicciones |
| 400-700 | Medium | **Recomendado** para análisis general |
| 700-900 | High | Alta confianza, interacciones bien validadas |
| 900-1000 | Highest | Máxima confianza, evidencia experimental |

#### Ejemplos de Uso

```bash
# 1. Análisis básico de un batch
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch01.faa \
  --out data/string/results \
  --species 4113

# 2. Alta confianza + enriquecimiento funcional
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch01.faa \
  --out data/string/results_high_conf \
  --min-score 700 \
  --get-enrichment

# 3. Procesar múltiples batches
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch*.faa \
  --out data/string/results_all \
  --min-score 400 \
  --get-enrichment

# 4. Red de interacciones físicas (solo)
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch01.faa \
  --out data/string/results_physical \
  --network-type physical \
  --min-score 700
```

#### Archivos de Salida

Para cada batch procesado se generan:

```
<OUTPUT_DIR>/
├── <batch_name>_network.tsv        # Red PPI (interacciones)
├── <batch_name>_enrichment.tsv     # Enriquecimiento funcional (si --get-enrichment)
└── <batch_name>_string_ids.tsv     # Mapeo de IDs originales → STRING IDs
```

**Formato de `*_network.tsv`**:
```
protein1    protein2    string_id1            string_id2            score
g5.t1       g12345.t1   4113.PGSC0003DM...   4113.PGSC0003DM...   785
```

**Formato de `*_enrichment.tsv`**:
```
category    term         description                    p_value    fdr
Process     GO:0006096   glycolytic process             1.2e-15    3.4e-13
KEGG        sot00010     Glycolysis / Gluconeogenesis   8.7e-10    2.1e-08
```

#### Limitaciones y Consideraciones

**⚠️ Rate Limiting**:
- El script implementa rate limiting automático (0.5s entre requests)
- Para análisis muy grandes (>10 batches), considera ejecutar en horarios de bajo tráfico
- Si recibes error 429 (Too Many Requests), el script esperará automáticamente

**⚠️ Tamaño de Red**:
- STRING limita networks a ~200 proteínas por query en la API
- Para batches más grandes, las interacciones se truncarán

**⚠️ Mapeo de IDs**:
- STRING usa identificadores propios (e.g., `4113.PGSC0003DMT400058594`)
- El archivo `*_string_ids.tsv` mantiene el mapeo a tus IDs originales

---

## Flujo de Trabajo Recomendado

### Fase 1: Exploración Inicial (Manual)

```bash
# 1. Generar archivos para web (ya hecho)
ls data/string/web_example/

# 2. Upload batch01 a STRING web interface
# - Ve a https://string-db.org/
# - Organism: Solanum tuberosum (4113)
# - Upload: data/string/web_example/creole_web_example_batch01.faa

# 3. Explorar:
# - Ajusta score de confianza visualmente
# - Identifica clusters/módulos funcionales
# - Exporta network como imagen para presentaciones
```

### Fase 2: Análisis Sistemático (API)

```bash
# 1. Procesar primer batch con enriquecimiento
python scripts/string/string_api_batch.py \
  --input data/string/web_example/creole_web_example_batch01.faa \
  --out data/string/results \
  --min-score 400 \
  --get-enrichment

# 2. Revisar resultados
head data/string/results/creole_web_example_batch01_network.tsv
head data/string/results/creole_web_example_batch01_enrichment.tsv

# 3. Si los resultados son buenos, escalar a más batches
# (Considera ejecutar batches individuales y consolidar después)
```

### Fase 3: Integración con Modelo GEM

```bash
# 1. Consolidar interacciones de todos los batches
cat data/string/results/*_network.tsv | grep -v "^protein1" | sort -t$'\t' -k5 -rn > data/string/all_interactions_sorted.tsv

# 2. Filtrar por score alto (≥700)
awk -F'\t' '$5 >= 700' data/string/all_interactions_sorted.tsv > data/string/high_confidence_interactions.tsv

# 3. Cruzar con GPRs del modelo actual
# (Script personalizado para tu pipeline)
python scripts/integrate_ppi_with_gem.py \
  --model models/current/creole_criolla_colombia_v1.3.xml \
  --ppi data/string/high_confidence_interactions.tsv \
  --out reports/ppi_gem_integration.txt
```

---

## Estadísticas de Datos Generados

### Batches para Web (Alta Confianza)
- **Ubicación**: `data/string/web_example/`
- **Filtros**: EC numbers requeridos, score ≥ 100
- **Total**: 11,832 proteínas
- **Batches**: 60 archivos × 200 proteínas
- **Cobertura**: 100% con EC numbers

### Batches Completos (Relevancia Metabólica)
- **Ubicación**: `data/string/batches/`
- **Filtros**: KEGG pathways/reactions/KO, score ≥ 50
- **Total**: 34,695 proteínas
- **Batches**: 70 archivos × 500 proteínas
- **Cobertura**:
  - 40.6% con EC numbers
  - 58.9% con KEGG pathways
  - 54.1% con GO terms

---

## Troubleshooting

### Error: "No proteins mapped to STRING IDs"

**Causa**: Los identificadores de proteínas no coinciden con la base de datos de STRING para *S. tuberosum*.

**Soluciones**:
1. Verifica que `--species 4113` es correcto
2. STRING puede no tener todos los genes de la variedad Phureja
3. Intenta con `--species 4113` (general potato) vs específico

### Error: "429 Too Many Requests"

**Causa**: Excediste el rate limit de la API de STRING.

**Solución**: El script esperará automáticamente 60s y reintentará. Para análisis grandes, considera:
- Ejecutar durante horarios de bajo tráfico (madrugada)
- Procesar batches individuales con pausas entre ejecuciones

### Warning: "Network truncated to 200 proteins"

**Causa**: STRING limita el tamaño de redes en queries API.

**Solución**:
- Reduce `--batch-size` a 100-150 al crear batches
- Ejecuta múltiples queries pequeñas y consolida resultados

---

## Recursos Adicionales

### Documentación de STRING
- **API Docs**: https://string-db.org/help/api/
- **FAQ**: https://string-db.org/help/faq/
- **Paper (2023)**: https://pubmed.ncbi.nlm.nih.gov/36370105/

### Datos de *Solanum tuberosum* en STRING
- **Taxonomy ID**: 4113
- **Organism Code**: `sot`
- **Total Proteins**: ~40,000 (depende de versión STRING)

### Integración con Otras Herramientas
- **Cytoscape**: Importa `*_network.tsv` para visualización avanzada
- **igraph/networkx**: Análisis de redes en Python/R
- **KEGG Mapper**: Mapea enriquecimiento a rutas metabólicas

---

## Contacto y Soporte

Para reportar issues o sugerencias sobre estos scripts:
- **Proyecto**: GEM_Creole - Genome-scale Metabolic Model para *S. tuberosum* Phureja
- **Ubicación**: `/home/oscar/Documents/Master/Thesis/GEM_Creole/`

---

## Changelog

### 2025-10-19 - Initial Release
- ✅ Implementado `filter_proteome_for_string.py`
- ✅ Implementado `string_api_batch.py`
- ✅ Generados batches de ejemplo para web y API
- ✅ Documentación completa en README
