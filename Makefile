SHELL := /usr/bin/env bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --no-builtin-rules

DATA := data
MAP  := $(DATA)/mappings
ANN  := $(DATA)/annotations
CACHE:= $(DATA)/cache
KC   := $(CACHE)/kegg
PIPE := $(CACHE)/pipeline/api_enrichment
SCRIPTS := scripts

# Entradas esperadas (producidas por tus AWK previos)
KO_EXP := $(MAP)/kegg_ko_expanded.csv
RXN_EXP:= $(MAP)/kegg_reactions_expanded.csv

# Salidas intermedias
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

kegg: ko-links rxn-details compounds

$(PIPE) $(KC)/rxn/flat $(KC)/compound/flat:
	mkdir -p $@

# -------- KOs
# KO list desde kegg_ko_expanded.csv (2a col), sin comillas
$(KO_LIST): $(KO_EXP) | $(PIPE)
	awk -F',' 'NR>0 {gsub(/"/,"",$2); if($2!="" && $2!="-") print $$2}' $< \
	| sort -u > $@

# KO -> (reaction/module/pathway) usando KEGG link
ko-links: $(KO_LIST) | $(PIPE)
	$(SCRIPTS)/pipeline/ko_links.sh $(KO_LIST) $(PIPE)

# -------- Reacciones
# RXN list = unificación de reacciones desde KO y tu CSV expandido (si existe)
$(RXN_LIST): $(RXN_FROM_KO) $(RXN_EXP) | $(PIPE)
	{ cut -f2 $(RXN_FROM_KO) 2>/dev/null | sed 's#^rn:##'; \
	  [ -f "$(RXN_EXP)" ] && awk -F',' '{gsub(/"/,"",$2); if(NR>0 && $$2!="") print $$2}' $(RXN_EXP) || true; } \
	| sort -u > $@

# Descarga detalles de reacciones en flat file (batch 10)
$(RXN_FLAT): $(RXN_LIST) | $(KC)/rxn/flat
	$(SCRIPTS)/pipeline/rxn_details.sh $(RXN_LIST) $@

# Parseo de flat -> CSVs (ecuación, ECs, compuestos, rclass)
rxn-details: $(RXN_FLAT)
	python3 $(SCRIPTS)/processing/parse_kegg_rxn_flat.py $(RXN_FLAT) \
	  --rxn-equation $(RXN_EQ) \
	  --rxn-ec $(RXN_EC) \
	  --rxn-compound-edges $(RXN_CEDG) \
	  --rclass-map $(RC_MAP)

# -------- Compounds
# Lista de compuestos Cxxxxx desde edges
$(CMP_LIST): $(RXN_CEDG) | $(PIPE)
	awk -F',' 'NR>1 {print $$3}' $(RXN_CEDG) | sort -u > $@

# Descarga compounds flat
$(CMP_FLAT): $(CMP_LIST) | $(KC)/compound/flat
	$(SCRIPTS)/pipeline/compounds_details.sh $(CMP_LIST) $@

# Parseo compounds flat -> CSV
compounds: $(CMP_FLAT)
	python3 $(SCRIPTS)/processing/parse_kegg_compound_flat.py $(CMP_FLAT) \
	  --out $(CMP_CSV)

clean:
	rm -f $(RXN_FLAT) $(CMP_FLAT) $(KO_LIST) $(RXN_LIST) $(CMP_LIST)
