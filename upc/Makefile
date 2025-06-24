# UPC FFT Makefile for Lab 204 Environment

# Default values
PTHREADS ?= 4
NODES ?= 4
THREADS = $(shell echo $$(($(PTHREADS) * $(NODES))))
INPUT ?= example1
SHARED_HEAP ?= 256M

# Executable names
TARGET_BUPC = fft_bupc
TARGET_GUPC = fft_gupc

# Source files
SOURCES = fft_upc.c

# Node file
NODEFILE = nodes

# Default target
all: setup compile-bupc

# Environment setup
setup:
	@echo "Setting up environment..."
	@echo "Make sure to run: source /opt/nfs/config/source_bupc.sh"
	@echo "Optionally: source /opt/nfs/config/source_mpich411.sh"

# Generate node file
nodes:
	/opt/nfs/config/station204_name_list.sh 1 16 > $(NODEFILE)
	@echo "Generated node file with available stations"

# Compilation with Berkeley UPC
compile-bupc: $(TARGET_BUPC)

$(TARGET_BUPC): $(SOURCES)
	upcc -bupc -network=udp -pthreads=$(PTHREADS) $(SOURCES) -o $(TARGET_BUPC) -lm

# Compilation with GNU UPC (alternative)
compile-gupc: $(TARGET_GUPC)

$(TARGET_GUPC): $(SOURCES)
	/opt/nfs/berkeley_upc-2021.4.0/bin/upcc -gupc -Wc,"-fPIE" -network=udp -pthreads=$(PTHREADS) $(SOURCES) -o $(TARGET_GUPC) -lm

# Debug compilation
compile-debug: fft_upc_debug

fft_upc_debug: fft_upc_debug.c
	upcc -bupc -network=udp -pthreads=$(PTHREADS) fft_upc_debug.c -o fft_upc_debug -lm

# Run debug version
run-debug: fft_upc_debug $(NODEFILE)
	@echo "Running debug version with $(THREADS) threads"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap $(SHARED_HEAP) -c $(PTHREADS) -N $(NODES) -n $(THREADS) ./fft_upc_debug $(INPUT)

# Run with Berkeley UPC
run: $(TARGET_BUPC) $(NODEFILE)
	@echo "Running with $(THREADS) threads ($(PTHREADS) threads/process × $(NODES) nodes)"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap $(SHARED_HEAP) -c $(PTHREADS) -N $(NODES) -n $(THREADS) ./$(TARGET_BUPC) $(INPUT)

# Run with GNU UPC
run-gupc: $(TARGET_GUPC) $(NODEFILE)
	@echo "Running with $(THREADS) threads ($(PTHREADS) threads/process × $(NODES) nodes)"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap $(SHARED_HEAP) -c $(PTHREADS) -N $(NODES) -n $(THREADS) ./$(TARGET_GUPC) $(INPUT)

# Generate example data
generate:
	python3 sygnalGenerate.py

plot:
	python3 fft_plot.py fft_output $(INPUT)


# Test with different configurations
test-small: $(TARGET_BUPC) $(NODEFILE)
	@echo "Testing with small configuration: 2 nodes, 2 threads/process"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap 128M -c 2 -N 2 -n 4 ./$(TARGET_BUPC) $(INPUT)

test-medium: $(TARGET_BUPC) $(NODEFILE)
	@echo "Testing with medium configuration: 4 nodes, 4 threads/process"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap 256M -c 4 -N 4 -n 16 ./$(TARGET_BUPC) $(INPUT)

test-large: $(TARGET_BUPC) $(NODEFILE)
	@echo "Testing with large configuration: 8 nodes, 4 threads/process"
	UPC_NODEFILE=$(NODEFILE) upcrun -shared-heap 512M -c 4 -N 8 -n 32 ./$(TARGET_BUPC) $(INPUT)

# Environment check
check-env:
	@echo "Checking UPC environment..."
	@which upcc > /dev/null 2>&1 && echo "✓ upcc compiler found" || echo "✗ upcc compiler not found - run setup first"
	@which upcrun > /dev/null 2>&1 && echo "✓ upcrun runtime found" || echo "✗ upcrun runtime not found"
	@test -f $(NODEFILE) && echo "✓ Node file exists" || echo "✗ Node file missing - run 'make nodes'"

# Clean compiled files
clean:
	rm -f $(TARGET_BUPC) $(TARGET_GUPC)

# Reset workspace (remove generated files too)
reset: clean
	rm -f example* fft_output plot.png $(NODEFILE)

# Help
help:
	@echo "UPC FFT Makefile for Lab 204 Environment"
	@echo ""
	@echo "Setup Commands:"
	@echo "  setup            - Show environment setup instructions"
	@echo "  nodes            - Generate node file for Lab 204"
	@echo "  check-env        - Check if UPC environment is properly set up"
	@echo ""
	@echo "Compilation:"
	@echo "  compile-bupc     - Compile with Berkeley UPC (default)"
	@echo "  compile-gupc     - Compile with GNU UPC"
	@echo "  compile-debug    - Compile debug version"
	@echo ""
	@echo "Execution:"
	@echo "  run              - Run with default parameters (Berkeley UPC)"
	@echo "  run-gupc         - Run with GNU UPC"
	@echo "  run-debug        - Run debug version"
	@echo "  run NODES=n      - Run with n nodes"
	@echo "  run PTHREADS=n   - Run with n threads per process"
	@echo "  run INPUT=file   - Run with custom input file"
	@echo ""
	@echo "Testing:"
	@echo "  test-small       - Test with 4 threads (2×2)"
	@echo "  test-medium      - Test with 16 threads (4×4)"
	@echo "  test-large       - Test with 32 threads (8×4)"
	@echo ""
	@echo "Utilities:"
	@echo "  generate         - Generate example input data"
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

.PHONY: all setup nodes compile-bupc compile-gupc compile-debug run run-gupc run-debug generate compare test-small test-medium test-large check-env clean reset help
