# Common part of Makefile for all reports

LATEXMK := latexmk
OUT_DIR := ../../output_files/reports/$(TARGET)
RC_FILE := ../common/.latexmkrc

.PHONY: FORCE_MAKE all clean cleanall

all: $(TARGET).pdf

%.pdf: %.tex FORCE_MAKE
	$(LATEXMK) -r $(RC_FILE) -output-directory=$(OUT_DIR) $<

clean:
	latexmk -c -output-directory=$(OUT_DIR)

cleanall:
	latexmk -C -output-directory=$(OUT_DIR)