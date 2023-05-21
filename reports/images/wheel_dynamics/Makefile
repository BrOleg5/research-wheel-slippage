TARGETS := lemma1.pdf lemma2.pdf

RM_COMMAND := 
ifeq ($(OS), Windows_NT)
	RM_COMMAND := del /Q
else
	RM_COMMAND := rm
endif

.PHONY: all clean

all: $(TARGETS)

%.pdf: %.svg
	inkscape --export-filename=$@ $<

clean:
	$(RM_COMMAND) *.pdf