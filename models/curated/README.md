# Curated Models - Modelos Curados

Este directorio contiene modelos que han pasado por procesos de curación y están listos para evaluación y uso en investigación.

## Modelos Disponibles

### Versión 1.0 - Modelo Final Optimizado
- `creole_v1.0_final.xml` - Modelo final tras 5 fases de mejora
- `creole_v1.0_final.consistency_report.json` - Reporte de consistencia estequiométrica  
- `creole_v1.0_final.consistency_summary.txt` - Resumen de análisis

### Modelos Previos
- `creole_with_gprs.xml` - Modelo con GPRs básicos
- `creole_with_gprs_annot.xml` - Versión con anotaciones adicionales

## Características del Modelo v1.0

- **1,085 reacciones** optimizadas
- **1,134 metabolitos** con anotaciones completas
- **2,117 genes** con anotaciones funcionales masivas (162K)
- **93.1% reacciones balanceadas** en masa
- **NGAM presente** y funcional
- **Biomasa optimizada** (GAM: 30→10 ATP)
- **84% consistencia estequiométrica**

## Validación
- ✅ FBA optimal (objetivo: 259.86)
- ✅ MEMOTE testing completo
- ✅ Balance masa/carga verificado
- ✅ Funcionalidad de biomasa validada
