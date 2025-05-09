CC = mpicc
CFLAGS = -Wall -O2 -std=c99
EXEC = fft_mpi_exec
EXEC_MPE = fft_mpe
SOURCE = fft.c
SOURCE_MPE = fft_mpe.c
INPUT = example
NUM_PROCESSES = 8 
NODES_SCRIPT = /opt/nfs/config/station204_name_list.sh
NODES_FILE = nodes
OUTPUT = fft_output
OUTPUT_PNG = plot.png
PYTHON_SCRIPT = fft_plot.py
PYTHON = python3
SHELL := /bin/bash

all: compile run

compile: $(EXEC)

compile_mpe: $(EXEC_MPE)

$(EXEC): $(SOURCE)
	$(CC) $(SOURCE) -o $(EXEC) $(CFLAGS) -lm

$(EXEC_MPE): $(SOURCE_MPE)
	$(CC) -g $(SOURCE_MPE) -o $(EXEC_MPE) $(CFLAGS) -I/opt/nfs/mpe2-2.4.9b/include -L/opt/nfs/mpe2-2.4.9b/lib -lmpe -lm

$(NODES_FILE):
	$(NODES_SCRIPT) 1 16 > $(NODES_FILE)

run: $(EXEC) $(INPUT) $(NODES_FILE)
	mpiexec -np $(NUM_PROCESSES) -f $(NODES_FILE) ./$(EXEC) $(INPUT)

run_mpe: $(EXEC_MPE) $(INPUT) $(NODES_FILE)
	mpiexec -np $(NUM_PROCESSES) -f $(NODES_FILE) ./$(EXEC_MPE) $(INPUT)

plot: run $(OUTPUT) $(INPUT)
	$(PYTHON) $(PYTHON_SCRIPT) $(OUTPUT) $(INPUT)

clean:
	rm -f $(EXEC) $(EXEC_MPE) *.o $(NODES_FILE) $(OUTPUT) $(OUTPUT_PNG) *.clog2 *.slog2 generated_*

convert:
	/opt/nfs/mpe2-2.4.9b/bin/clog2TOslog2 mpe_log.clog2

jumpshot:
	/opt/nfs/mpe2-2.4.9b/bin/jumpshot mpe_log.slog2

reset: clean
	git checkout . # Przywraca zmiany w śledzonych plikach
	git clean -fd # Usuwa nieśledzone pliki i katalogi (ostrożnie z tą komendą!)

.PHONY: all compile run clean reset
