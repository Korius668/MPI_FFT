# UPC FFT Makefile

# Default values
PTHREADS ?= 4
NODES ?= 4
THREADS = $(shell echo $$(($(PTHREADS) * $(NODES))))
INPUT ?= example
SHARED_HEAP ?= 256M

# Executable names
TARGET_BUPC = fft_bupc

# Source files
SOURCES = fft_upc.c

# Node file
NODEFILE = nodes

# Default target
all: setup compile run

# Environment setup
setup:
	@echo "Setting up environment..."
	@echo "Make sure to run: source /opt/nfs/config/source_bupc.sh"

# Generate node file
nodes:
	/opt/nfs/config/station204_name_list.sh 1 16 > $(NODEFILE)
	@echo "Generated node file with available stations"

# Compilation with Berkeley UPC
compile: $(TARGET_BUPC)

$(TARGET_BUPC): $(SOURCES)
	upcc -bupc -network=udp -pthreads=$(PTHREADS) $(SOURCES) -o $(TARGET_BUPC) -lm

# Run with Berkeley UPC
run: $(TARGET_BUPC) $(NODEFILE)
	@echo "Running with $(THREADS) threads ($(PTHREADS) threads/process × $(NODES) nodes)"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap $(SHARED_HEAP) -c $(PTHREADS) -N $(NODES) -n $(THREADS) ./$(TARGET_BUPC) $(INPUT)

plot:
	python3 fft_plot.py fft_output $(INPUT)

# Clean compiled files
clean:
	rm -f $(TARGET_BUPC) $(TARGET_GUPC)

# Reset workspace (remove generated files too)
reset: clean
	rm -f plot.png $(NODEFILE) fft_output generate_*

# Help
help:
	@echo "UPC FFT Makefile"
	@echo ""
	@echo "Setup Commands:"
	@echo "  setup            - Show environment setup instructions"
	@echo "  nodes            - Generate node file for Lab 204"
	@echo ""
	@echo "Compilation:"
	@echo "  compile    - Compile with Berkeley UPC (default)"
	@echo ""
	@echo "Execution:"
	@echo "  run              - Run with default parameters (Berkeley UPC)"
	@echo "  run NODES=n      - Run with n nodes"
	@echo "  run PTHREADS=n   - Run with n threads per process"
	@echo "  run INPUT=file   - Run with custom input file"
	@echo ""
	@echo "Utilities:"
	@echo "  compare          - Compare results with expected output"
	@echo "  clean            - Remove compiled files"
	@echo "  reset            - Remove all generated files"
	@echo ""
	@echo "Parameters:"
	@echo "  PTHREADS=$(PTHREADS)     - Threads per process"
	@echo "  NODES=$(NODES)        - Number of nodes"
	@echo "  THREADS=$(THREADS)       - Total threads (PTHREADS × NODES)"
	@echo "  INPUT=$(INPUT)      - Input file"
	@echo "  SHARED_HEAP=$(SHARED_HEAP)  - Shared heap size"

.PHONY: all setup nodes compile run compare clean reset help
