#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <upc.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef OUTPUT
#define OUTPUT "fft_output"
#endif

// Define our own complex number structure for UPC compatibility
typedef struct {
    double real;
    double imag;
} complex_t;

// Complex number operations
complex_t complex_add(complex_t a, complex_t b) {
    complex_t result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

complex_t complex_sub(complex_t a, complex_t b) {
    complex_t result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

complex_t complex_mul(complex_t a, complex_t b) {
    complex_t result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

complex_t complex_exp(double angle) {
    complex_t result;
    result.real = cos(angle);
    result.imag = sin(angle);
    return result;
}

complex_t make_complex(double real, double imag) {
    complex_t result;
    result.real = real;
    result.imag = imag;
    return result;
}

// Shared variables
shared int N;
shared int Fs;

int main(int argc, char *argv[]) {
    printf("Thread %d of %d: Starting debug version\n", MYTHREAD, THREADS);
    
    if (MYTHREAD == 0) {
        if (argc != 2) {
            fprintf(stderr, "Usage: %s <signal_file>\n", argv[0]);
            exit(1);
        }
        
        FILE *infile = fopen(argv[1], "r");
        if (infile == NULL) {
            perror("Cannot open input file");
            exit(1);
        }
        
        int temp_Fs;
        if (fscanf(infile, "%d", &temp_Fs) != 1) {
            fprintf(stderr, "Error reading sampling frequency\n");
            fclose(infile);
            exit(1);
        }
        Fs = temp_Fs;
        printf("Thread 0: Read Fs = %d\n", Fs);
        
        // Count samples
        double val;
        int count = 0;
        while (fscanf(infile, "%lf", &val) == 1) {
            count++;
        }
        fclose(infile);
        
        N = count;
        printf("Thread 0: Read %d samples\n", N);
        
        // Check divisibility
        if (N % THREADS != 0) {
            printf("Thread 0: Padding N from %d to %d\n", N, ((N + THREADS - 1) / THREADS) * THREADS);
            N = ((N + THREADS - 1) / THREADS) * THREADS;
        }
    }
    
    upc_barrier;
    
    printf("Thread %d: N = %d, local_N = %d\n", MYTHREAD, N, N/THREADS);
    
    // Test shared memory allocation
    shared [*] complex_t *test_array = (shared [*] complex_t *) upc_all_alloc(THREADS, (N/THREADS) * sizeof(complex_t));
    if (test_array == NULL) {
        printf("Thread %d: Shared memory allocation FAILED\n", MYTHREAD);
        exit(1);
    } else {
        printf("Thread %d: Shared memory allocation SUCCESS\n", MYTHREAD);
    }
    
    upc_barrier;
    
    // Test writing to shared memory
    int local_N = N / THREADS;
    for (int i = 0; i < local_N; i++) {
        test_array[MYTHREAD * local_N + i] = make_complex(MYTHREAD, i);
    }
    
    upc_barrier;
    
    // Test reading from shared memory
    if (MYTHREAD == 0) {
        printf("Thread 0: Testing shared memory access:\n");
        for (int t = 0; t < THREADS && t < 4; t++) {
            for (int i = 0; i < local_N && i < 4; i++) {
                complex_t val = test_array[t * local_N + i];
                printf("  [%d][%d] = %.1f + %.1fi\n", t, i, val.real, val.imag);
            }
        }
    }
    
    upc_free(test_array);
    
    printf("Thread %d: Debug test completed successfully\n", MYTHREAD);
    
    return 0;
}
