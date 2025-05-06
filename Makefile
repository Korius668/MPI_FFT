CC = mpicc
CFLAGS = -Wall -O2 -std=c99
EXEC = fft_mpi_exec
SOURCE = fft.c
INPUT = sygnal
NUM_PROCESSES = 8 
NODES_SCRIPT = /opt/nfs/config/station204_name_list.sh
NODES_FILE = nodes
OUTPUT = fft_output
OUTPUT_PNG = fft_spectrum.png
PYTHON_SCRIPT = fft_plot.py
PYTHON = python3

all: compile run

compile: $(EXEC)

$(EXEC): $(SOURCE)
	$(CC) $(SOURCE) -o $(EXEC) $(CFLAGS) -lm

$(NODES_FILE):
	$(NODES_SCRIPT) 1 16 > $(NODES_FILE)

run: $(EXEC) $(INPUT) $(NODES_FILE)
	mpiexec -np $(NUM_PROCESSES) -f $(NODES_FILE) ./$(EXEC) $(INPUT)

plot: run $(OUTPUT) $(INPUT)
	$(PYTHON) $(PYTHON_SCRIPT) $(OUTPUT) $(INPUT)

clean:
	rm -f $(EXEC) *.o $(NODES_FILE) $(OUTPUT) $(OUTPUT_PNG)

reset: clean
	git checkout . # Przywraca zmiany w śledzonych plikach
	git clean -fd # Usuwa nieśledzone pliki i katalogi (ostrożnie z tą komendą!)

.PHONY: all compile run clean reset
