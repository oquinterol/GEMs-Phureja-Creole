SHELL := /bin/bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --no-builtin-rules

DATA := data
MAP  := $(DATA)/mappings
ANN  := $(DATA)/annotations
CACHE:= $(DATA)/cache
KC   := $(CACHE)/kegg
PIPE := $(CACHE)/pipeline/api_enrichment
SCRIPTS := scripts

# Entradas (formato Opción A: una fila por id, col 2 con lista separada por comas)
KO_GRP  := $(MAP)/kegg_ko.csv
RXN_GRP := $(MAP)/kegg_reactions.csv

# Intermedios
KO_LIST := $(PIPE)/ko.list
RXN_FROM_KO := $(PIPE)/ko2rxn.tsv
MOD_FROM_KO := $(PIPE)/ko2mod.tsv
MAP_FROM_KO := $(PIPE)/ko2map.tsv
RXN_LIST := $(PIPE)/rxn.list
RXN_FLAT := $(KC)/rxn/flat/rxn_flat.txt

# Salidas finales (CSV)
RXN_EQ   := $(MAP)/rxn_equation.csv
RXN_EC   := $(MAP)/rxn_ec.csv
RXN_CEDG := $(MAP)/rxn_compound_edges.csv
RC_MAP   := $(MAP)/rclass_map.csv
CMP_LIST := $(PIPE)/compounds.list
CMP_FLAT := $(KC)/compound/flat/compound_flat.txt
CMP_CSV  := $(MAP)/compounds.csv

.PHONY: all kegg ko-links rxn-details compounds clean

all: kegg

# “kegg” como alias que exige todos los productos finales (idempotente)
kegg: $(RXN_EQ) $(RXN_EC) $(RXN_CEDG) $(RC_MAP) $(CMP_CSV)

# Alias de conveniencia
rxn-details: $(RXN_EQ) $(RXN_EC) $(RXN_CEDG) $(RC_MAP)
compounds:   $(CMP_CSV)
ko-links:    $(RXN_FROM_KO) $(MOD_FROM_KO) $(MAP_FROM_KO)

$(PIPE) $(KC)/rxn/flat $(KC)/compound/flat:
	mkdir -p $@

# -------- KO list desde kegg_ko.csv (col 2 con lista separada por comas, Opción A)
$(KO_LIST): $(KO_GRP) | $(PIPE)
	awk -F',' 'NR>1 { \
	  s=$$2; gsub(/^"|"$$/,"",s); \
	  n=split(s,a,/,/); \
	  for(i=1;i<=n;i++){ \
	    gsub(/^ +| +$$/,"",a[i]); \
	    gsub(/^ko:/,"",a[i]);            # <--- quita prefijo ko: si viene en el CSV
	    if(a[i]!="" && a[i]!="-") print a[i]; \
	  } \
	}' $< | sort -u > $@

# -------- KO -> (reaction/module/pathway) via KEGG link (targets reales)
$(RXN_FROM_KO) $(MOD_FROM_KO) $(MAP_FROM_KO): $(KO_LIST) | $(PIPE)
	@echo ">>> Ejecutando ko_links.sh..."
	$(SCRIPTS)/pipeline/ko_links.sh $(KO_LIST) $(PIPE)
	@# Asegura timestamps de los tres outputs
	@touch $(RXN_FROM_KO) $(MOD_FROM_KO) $(MAP_FROM_KO)
	@echo ">>> ko_links.sh finalizado: $(RXN_FROM_KO) $(MOD_FROM_KO) $(MAP_FROM_KO)"

# -------- RXN list = ko2rxn.tsv + (opcional) kegg_reactions.csv (Opción A)
$(RXN_LIST): $(RXN_FROM_KO) $(RXN_GRP) | $(PIPE)
	{ \
	  cut -f2 $(RXN_FROM_KO) 2>/dev/null | sed 's#^rn:##'; \
	  awk -F',' 'NR>1 { \
	    s=$$2; gsub(/^"|"$$/,"",s); \
	    n=split(s,a,/,/); \
	    for(i=1;i<=n;i++){ gsub(/^ +| +$$/,"",a[i]); if(a[i]!="") print a[i]; } \
	  }' $(RXN_GRP) 2>/dev/null || true; \
	} | sort -u > $@

# -------- Descargar detalles de reacciones (flat)
$(RXN_FLAT): $(RXN_LIST) | $(KC)/rxn/flat
	@echo ">>> Descargando detalles de reacciones (flat)..."
	$(SCRIPTS)/pipeline/rxn_details.sh $(RXN_LIST) $@
	@echo ">>> rxn_flat listo: $@"

# -------- Parsear flat -> CSVs (ecuación, ECs, compuestos, rclass) -> targets reales
$(RXN_EQ) $(RXN_EC) $(RXN_CEDG) $(RC_MAP): $(RXN_FLAT)
	@echo ">>> Ejecutando parse_kegg_rxn_flat.py sobre $(RXN_FLAT)..."
	python3 -u $(SCRIPTS)/processing/parse_kegg_rxn_flat.py $(RXN_FLAT) \
	  --rxn-equation $(RXN_EQ) \
	  --rxn-ec $(RXN_EC) \
	  --rxn-compound-edges $(RXN_CEDG) \
	  --rclass-map $(RC_MAP)
	@echo ">>> parse_kegg_rxn_flat.py finalizado:"
	@ls -lh $(RXN_EQ) $(RXN_EC) $(RXN_CEDG) $(RC_MAP) || true

# -------- Compounds: lista desde edges
$(CMP_LIST): $(RXN_CEDG) | $(PIPE)
	awk -F',' 'NR>1 {print $$3}' $(RXN_CEDG) | sort -u > $@

# -------- Descargar compounds flat
$(CMP_FLAT): $(CMP_LIST) | $(KC)/compound/flat
	@echo ">>> Descargando compounds (flat)..."
	$(SCRIPTS)/pipeline/compounds_details.sh $(CMP_LIST) $@
	@echo ">>> compound_flat listo: $@"

# -------- Parsear compounds flat -> CSV (target real)
$(CMP_CSV): $(CMP_FLAT)
	@echo ">>> Ejecutando parse_kegg_compound_flat.py..."
	python3 -u $(SCRIPTS)/processing/parse_kegg_compound_flat.py $(CMP_FLAT) \
	  --out $(CMP_CSV)
	@echo ">>> parse_kegg_compound_flat.py finalizado:"
	@ls -lh $(CMP_CSV) || true

clean:
	rm -f $(RXN_FLAT) $(CMP_FLAT) $(KO_LIST) $(RXN_LIST) $(CMP_LIST)
