# Current Model - Versión Oficial Actual

Este directorio contiene únicamente la versión oficial actual del modelo metabólico de *Solanum tuberosum* cultivar 'Criolla Colombia'.

## Modelo Actual

### creole_criolla_colombia_v1.0.xml
**Modelo Metabólico Oficial - Versión 1.0**

- **Organismo**: *Solanum tuberosum* cultivar 'Criolla Colombia'  
- **Nivel SBML**: Level 3 Version 2
- **Fecha**: Septiembre 2025
- **Estado**: Validado y optimizado

#### Especificaciones Técnicas:
- **1,085 reacciones** metabólicas optimizadas
- **1,134 metabolitos** con anotaciones completas  
- **2,117 genes** con anotaciones funcionales masivas
- **93.1% reacciones balanceadas** en masa y carga
- **NGAM**: Presente (ATP maintenance)
- **Biomasa**: Optimizada (GAM = 10 ATP)
- **Consistencia**: 84% estequiométrica

#### Validación:
- ✅ **FBA Optimal**: Objetivo 259.86
- ✅ **MEMOTE Testing**: Mejoras significativas en todas las categorías
- ✅ **Balance**: Masa y carga verificados
- ✅ **Funcionalidad**: Biomasa validada experimentalmente

#### Anotaciones Funcionales:
- **162,583 anotaciones** de eggNOG-mapper integradas
- **386 vías KEGG** organizadas
- **GO terms**, **EC numbers**, **KEGG pathways**
- **2,053 genes** enriquecidos funcionalmente

## Archivos de Soporte

- `creole_criolla_colombia_v1.0.consistency_report.json` - Análisis de consistencia detallado
- `creole_criolla_colombia_v1.0.consistency_summary.txt` - Resumen ejecutivo

## Uso Recomendado

Este modelo está listo para:
- Análisis de balance de flujo (FBA)
- Estudios de ingeniería metabólica  
- Análisis de redes metabólicas
- Investigación en biotecnología de papa
- Estudios de sistemas biológicos

## Historial

Para versiones anteriores y desarrollo, consultar:
- `models/drafts/` - Versiones de trabajo
- `models/curated/` - Versiones curadas

## Citación

Modelo desarrollado mediante proceso iterativo de 5 fases con integración masiva de anotaciones funcionales y optimización sistemática según estándares MEMOTE.
