# ==== Paths ====
RAW_DIR       := data/raw
PROC_DIR      := data/processed
FIG_DIR       := results/figs
TAB_DIR       := results/tables
REPORT_DIR    := report
SCRIPTS_DIR   := scripts

# ==== Default target ====
.PHONY: all
all: help

# ==== Setup (install R packages) ====
.PHONY: setup
setup:
	Rscript $(SCRIPTS_DIR)/00_setup.R

# ==== Task 1: hypercube + K-means + clusGap ====
T1_DATA  := $(PROC_DIR)/t1_sim.rds
T1_TAB   := $(TAB_DIR)/t1_gap_summary.csv
T1_FIG   := $(FIG_DIR)/t1_gap_curve.png

$(PROC_DIR)/t1_sim.rds: $(SCRIPTS_DIR)/task1_generate.R
	mkdir -p $(PROC_DIR)
	Rscript $(SCRIPTS_DIR)/task1_generate.R --out $(PROC_DIR)/t1_sim.rds

$(TAB_DIR)/t1_gap_summary.csv $(FIG_DIR)/t1_gap_curve.png: $(SCRIPTS_DIR)/task1_gap_kmeans.R $(PROC_DIR)/t1_sim.rds $(SCRIPTS_DIR)/utils_gap_helpers.R
	mkdir -p $(TAB_DIR) $(FIG_DIR)
	Rscript $(SCRIPTS_DIR)/task1_gap_kmeans.R \
	  --in $(PROC_DIR)/t1_sim.rds \
	  --tab $(TAB_DIR)/t1_gap_summary.csv \
	  --fig $(FIG_DIR)/t1_gap_curve.png \
	  --kmax_cap 12 \
	  --B 50

.PHONY: t1
t1: $(T1_TAB) $(T1_FIG)

# ==== Help ====
.PHONY: help
help:
	@echo "make setup    # Install R packages"
	@echo "make t1       # Run Task 1: simulate + clusGap + plot"
	@echo "make          # Show this help"
