#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <upc.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef OUTPUT
#define OUTPUT "fft_output"
#endif

typedef struct {
    double real;
    double imag;
} complex_t;

complex_t complex_add(complex_t a, complex_t b) {
    complex_t result = {a.real + b.real, a.imag + b.imag};
    return result;
}

complex_t complex_sub(complex_t a, complex_t b) {
    complex_t result = {a.real - b.real, a.imag - b.imag};
    return result;
}

complex_t complex_mul(complex_t a, complex_t b) {
    complex_t result = {
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    };
    return result;
}

complex_t complex_exp(double angle) {
    complex_t result = {cos(angle), sin(angle)};
    return result;
}

complex_t make_complex(double real, double imag) {
    complex_t result = {real, imag};
    return result;
}

shared int N;
shared int Fs;
shared [] complex_t *shared_signal;

int is_power_of_two(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

complex_t* pad_to_power_of_two(complex_t *signal, int n, int *padded_n) {
    if (is_power_of_two(n)) {
        *padded_n = n;
        complex_t *padded_signal = malloc(n * sizeof(complex_t));
        for (int i = 0; i < n; i++) padded_signal[i] = signal[i];
        return padded_signal;
    } else {
        int next_power = 1;
        while (next_power < n) next_power <<= 1;
        *padded_n = next_power;
        complex_t *padded_signal = malloc(next_power * sizeof(complex_t));
        for (int i = 0; i < n; i++) padded_signal[i] = signal[i];
        for (int i = n; i < next_power; i++) padded_signal[i] = make_complex(0.0, 0.0);
        return padded_signal;
    }
}

void swap(complex_t *a, complex_t *b) {
    complex_t tmp = *a; *a = *b; *b = tmp;
}

void bit_reverse_permutation(complex_t *x, int N) {
    int bits = 0;
    while ((1 << bits) < N) bits++;
    for (int i = 0; i < N; i++) {
        int rev = 0, tmp = i;
        for (int j = 0; j < bits; j++) {
            rev = (rev << 1) | (tmp & 1);
            tmp >>= 1;
        }
        if (rev > i) swap(&x[i], &x[rev]);
    }
}

void fft_compute_stages_upc(int N, int local_N) {
    int stages = 0;
    while ((1 << stages) < N) stages++;
    for (int s = 1; s <= stages; s++) {
        int m = 1 << s;
        int m2 = m >> 1;
        complex_t Wm = complex_exp(-2.0 * M_PI / m);
        for (int k = 0; k < N; k += m) {
            complex_t W = make_complex(1.0, 0.0);
            for (int j = 0; j < m2; j++) {
                int idx1 = k + j;
                int idx2 = k + j + m2;
                if (upc_threadof(&shared_signal[idx1]) == MYTHREAD ||
                    upc_threadof(&shared_signal[idx2]) == MYTHREAD) {
                    complex_t u = shared_signal[idx1];
                    complex_t t = complex_mul(W, shared_signal[idx2]);
                    if (upc_threadof(&shared_signal[idx1]) == MYTHREAD)
                        shared_signal[idx1] = complex_add(u, t);
                    if (upc_threadof(&shared_signal[idx2]) == MYTHREAD)
                        shared_signal[idx2] = complex_sub(u, t);
                }
                W = complex_mul(W, Wm);
            }
        }
        upc_barrier;
    }
}

int main(int argc, char *argv[]) {
    int local_N;
    complex_t *signal_full = NULL;
    complex_t *fft_result = NULL;
    struct timeval tv_start, tv_end;
    double start_time, end_time;

    if (MYTHREAD == 0) {
        if (argc != 2) {
            fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
            exit(1);
        }
        FILE *f = fopen(argv[1], "r");
        if (!f) { perror("open"); exit(1); }

        int Fs_tmp, count = 0, cap = 1024;
        fscanf(f, "%d", &Fs_tmp);
        Fs = Fs_tmp;
        signal_full = malloc(cap * sizeof(complex_t));
        double val;
        while (fscanf(f, "%lf", &val) == 1) {
            if (count >= cap) {
                cap *= 2;
                signal_full = realloc(signal_full, cap * sizeof(complex_t));
            }
            signal_full[count++] = make_complex(val, 0.0);
        }
        fclose(f);

        int padded_N;
        complex_t *padded = pad_to_power_of_two(signal_full, count, &padded_N);
        if (padded != signal_full) free(signal_full);
        signal_full = padded;

        if (padded_N % THREADS != 0) {
            fprintf(stderr, "Error: Padded length must be divisible by thread count\n");
            exit(1);
        }
        N = padded_N;
        bit_reverse_permutation(signal_full, N);
        fft_result = malloc(N * sizeof(complex_t));
    }

    upc_barrier;
    local_N = N / THREADS;

    // Dynamic shared memory allocation after N is known
    shared_signal = (shared [] complex_t *) upc_alloc(N * sizeof(complex_t));
    if (shared_signal == NULL) {
        fprintf(stderr, "Thread %d: shared_signal allocation failed\n", MYTHREAD);
        upc_global_exit(1);
    }

    // Data distribution
    if (MYTHREAD == 0) {
        for (int i = 0; i < N; i++) shared_signal[i] = signal_full[i];
        free(signal_full);
    }
    upc_barrier;

    if (MYTHREAD == 0) gettimeofday(&tv_start, NULL);
    fft_compute_stages_upc(N, local_N);
    if (MYTHREAD == 0) gettimeofday(&tv_end, NULL);

    // Result collection
    if (MYTHREAD == 0) {
        for (int i = 0; i < N; i++) fft_result[i] = shared_signal[i];
        start_time = tv_start.tv_sec + tv_start.tv_usec / 1e6;
        end_time = tv_end.tv_sec + tv_end.tv_usec / 1e6;
        printf("Execution time: %f seconds\n", end_time - start_time);

        FILE *out = fopen(OUTPUT, "w");
        fprintf(out, "%d %d\n", N, Fs);
        for (int i = 0; i < N; i++)
            fprintf(out, "%.8f %.8f\n", fft_result[i].real, fft_result[i].imag);
        fclose(out);
        free(fft_result);
    }

    upc_free(shared_signal);
    return 0;
}


// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <time.h>
// #include <sys/time.h>
// #include <upc.h>

// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

// #ifndef OUTPUT
// #define OUTPUT "fft_output"
// #endif

// // Define our own complex number structure for UPC compatibility
// typedef struct {
//     double real;
//     double imag;
// } complex_t;

// // Complex number operations
// complex_t complex_add(complex_t a, complex_t b) {
//     complex_t result;
//     result.real = a.real + b.real;
//     result.imag = a.imag + b.imag;
//     return result;
// }

// complex_t complex_sub(complex_t a, complex_t b) {
//     complex_t result;
//     result.real = a.real - b.real;
//     result.imag = a.imag - b.imag;
//     return result;
// }

// complex_t complex_mul(complex_t a, complex_t b) {
//     complex_t result;
//     result.real = a.real * b.real - a.imag * b.imag;
//     result.imag = a.real * b.imag + a.imag * b.real;
//     return result;
// }

// complex_t complex_exp(double angle) {
//     complex_t result;
//     result.real = cos(angle);
//     result.imag = sin(angle);
//     return result;
// }

// complex_t make_complex(double real, double imag) {
//     complex_t result;
//     result.real = real;
//     result.imag = imag;
//     return result;
// }

// // Shared variables - properly declared
// shared int N;
// shared int Fs;
// shared complex_t *shared_signal;

// int is_power_of_two(int n) {
//     return (n > 0) && ((n & (n - 1)) == 0);
// }

// // Function to pad signal with zeros to nearest power of 2
// complex_t* pad_to_power_of_two(complex_t *signal, int n, int *padded_n) {
//     if (is_power_of_two(n)) {
//         *padded_n = n;
//         complex_t *padded_signal = (complex_t*)malloc(n * sizeof(complex_t));
//         if (padded_signal == NULL) {
//             fprintf(stderr, "Thread %d: Memory allocation error\n", MYTHREAD);
//             exit(1);
//         }
//         for (int i = 0; i < n; i++) {
//             padded_signal[i] = signal[i];
//         }
//         return padded_signal;
//     } else {
//         int next_power_of_two = 1;
//         while (next_power_of_two < n) {
//             next_power_of_two <<= 1;
//         }
//         *padded_n = next_power_of_two;
//         complex_t *padded_signal = (complex_t*)malloc(next_power_of_two * sizeof(complex_t));
//         if (padded_signal == NULL) {
//             fprintf(stderr, "Thread %d: Memory allocation error\n", MYTHREAD);
//             exit(1);
//         }
//         for (int i = 0; i < n; i++) {
//             padded_signal[i] = signal[i];
//         }
//         for (int i = n; i < next_power_of_two; i++) {
//             padded_signal[i] = make_complex(0.0, 0.0);
//         }
//         return padded_signal;
//     }
// }

// // Helper function to swap two complex numbers
// void swap(complex_t *a, complex_t *b) {
//     complex_t temp = *a;
//     *a = *b;
//     *b = temp;
// }

// // Bit-reversal permutation (performed by thread 0)
// void bit_reverse_permutation(complex_t *x, int N) {
//     int bits = 0;
//     while ((1 << bits) < N) {
//         bits++;
//     }

//     for (int i = 0; i < N; i++) {
//         int rev = 0;
//         int temp = i;
//         for (int j = 0; j < bits; j++) {
//             rev = (rev << 1) | (temp & 1);
//             temp >>= 1;
//         }
//         if (rev > i) {
//             swap(&x[i], &x[rev]);
//         }
//     }
// }

// // UPC FFT computation using shared memory
// void fft_compute_stages_upc(int N, int local_N) {
//     int stages = 0;
//     while ((1 << stages) < N) {
//         stages++;
//     }

//     complex_t W, Wm, t, u;
//     int m, m2;

//     // Iterate through FFT stages
//     for (int s = 1; s <= stages; s++) {
//         m = 1 << s;  // Size of current sub-FFT
//         m2 = m >> 1; // Half size
//         Wm = complex_exp(-2.0 * M_PI / m); // Twiddle factor for stage

//         // Iterate through butterflies in stage
//         for (int k = 0; k < N; k += m) {
//             W = make_complex(1.0, 0.0); // Reset twiddle factor for each butterfly group
//             for (int j = 0; j < m2; j++) {
//                 int idx1 = k + j;
//                 int idx2 = k + j + m2;

//                 // Determine which threads own the needed indices
//                 int thread1 = idx1 / local_N;
//                 int thread2 = idx2 / local_N;
//                 int local_idx1 = idx1 % local_N;
//                 int local_idx2 = idx2 % local_N;

//                 if (MYTHREAD == thread1 && MYTHREAD == thread2) {
//                     // Both indices are local - standard butterfly
//                     t = complex_mul(W, shared_signal[MYTHREAD * local_N + local_idx2]);
//                     u = shared_signal[MYTHREAD * local_N + local_idx1];
//                     shared_signal[MYTHREAD * local_N + local_idx1] = complex_add(u, t);
//                     shared_signal[MYTHREAD * local_N + local_idx2] = complex_sub(u, t);
//                 } else if (MYTHREAD == thread1) {
//                     // This thread has idx1, needs idx2 from thread2
//                     complex_t val_idx2 = shared_signal[thread2 * local_N + local_idx2];
                    
//                     // Compute butterfly using remote value
//                     t = complex_mul(W, val_idx2);
//                     u = shared_signal[MYTHREAD * local_N + local_idx1];
//                     shared_signal[MYTHREAD * local_N + local_idx1] = complex_add(u, t);
//                 } else if (MYTHREAD == thread2) {
//                     // This thread has idx2, needs idx1 from thread1
//                     complex_t val_idx1 = shared_signal[thread1 * local_N + local_idx1];
                    
//                     // Compute butterfly using remote value
//                     t = complex_mul(W, shared_signal[MYTHREAD * local_N + local_idx2]);
//                     u = val_idx1;
//                     shared_signal[MYTHREAD * local_N + local_idx2] = complex_sub(u, t);
//                 }
                
//                 // Update twiddle factor for next butterfly in group
//                 W = complex_mul(W, Wm);
//             }
//         }
        
//         // Synchronization barrier after each stage
//         upc_barrier;
//     }
// }

// int main(int argc, char *argv[]) {
//     int local_N; // Local signal size per thread
//     complex_t *signal_full = NULL; // Full signal (only on thread 0)
//     complex_t *fft_result = NULL;   // FFT result (only on thread 0)
//     char *filename = NULL;
//     FILE *infile = NULL;
//     struct timeval tv_start, tv_end;
//     double start_time, end_time;
//     int temp_N = 0;
    
//     if (MYTHREAD == 0) {
//         printf("Thread 0: Starting UPC FFT with %d threads\n", THREADS);
//     }
    
//     // Step 1: Data reading and distribution (Thread 0)
//     if (MYTHREAD == 0) {
//         if (argc != 2) {
//             fprintf(stderr, "Usage: %s <signal_file>\n", argv[0]);
//             exit(1);
//         }
//         filename = argv[1];
//         infile = fopen(filename, "r");
//         if (infile == NULL) {
//             perror("Cannot open input file");
//             exit(1);
//         }

//         int temp_Fs;
//         if (fscanf(infile, "%d", &temp_Fs) != 1) {
//             fprintf(stderr, "Error reading sampling frequency from file.\n");
//             fclose(infile);
//             exit(1);
//         }
//         Fs = temp_Fs;
//         printf("Thread 0: Read sampling frequency: %d Hz.\n", Fs);

//         // Read data and count N
//         double real_val;
//         int capacity = 1024;
//         signal_full = (complex_t *)malloc(capacity * sizeof(complex_t));
//         if (!signal_full) {
//             perror("Memory allocation error (signal_full)");
//             fclose(infile);
//             exit(1);
//         }

//         while (fscanf(infile, "%lf", &real_val) == 1) {
//             if (temp_N >= capacity) {
//                 capacity *= 2;
//                 complex_t *temp = (complex_t *)realloc(signal_full, capacity * sizeof(complex_t));
//                 if (!temp) {
//                     perror("Memory reallocation error (signal_full)");
//                     free(signal_full);
//                     fclose(infile);
//                     exit(1);
//                 }
//                 signal_full = temp;
//             }
//             signal_full[temp_N] = make_complex(real_val, 0.0);
//             temp_N++;
//         }
//         fclose(infile);

//         // Check N conditions
//         if (temp_N == 0) {
//             fprintf(stderr, "Error: Input file is empty.\n");
//             free(signal_full);
//             exit(1);
//         }
        
//         // Check if N is power of 2
//         int padded_N;
//         complex_t *padded_signal_full = pad_to_power_of_two(signal_full, temp_N, &padded_N);
//         if (padded_signal_full != signal_full) {
//             free(signal_full);
//             signal_full = padded_signal_full;
//             temp_N = padded_N;
//             printf("Thread 0: Signal padded with zeros to length %d (nearest power of 2).\n", temp_N);
//         } else {
//             printf("Thread 0: Number of samples (%d) is already a power of 2.\n", temp_N);
//         }

//         // Check if N is divisible by number of threads
//         if (temp_N % THREADS != 0) {
//             fprintf(stderr, "Error: Number of samples (%d) must be divisible by number of threads (%d).\n", temp_N, THREADS);
//             free(signal_full);
//             exit(1);
//         }

//         N = temp_N;
//         printf("Thread 0: Read %d samples from file '%s'. Number of threads: %d.\n", N, filename, THREADS);

//         // Perform bit-reversal permutation on full signal before distribution
//         bit_reverse_permutation(signal_full, N);

//         // Allocate memory for final result
//         fft_result = (complex_t *)malloc(N * sizeof(complex_t));
//         if (!fft_result) {
//             perror("Memory allocation error (fft_result)");
//             free(signal_full);
//             exit(1);
//         }
//     }

//     // Synchronize all threads before proceeding
//     upc_barrier;

//     // All threads now know N
//     local_N = N / THREADS;
    
//     if (MYTHREAD == 0) {
//         printf("Thread 0: Local N per thread: %d\n", local_N);
//     }

//     // Step 2: Allocate shared memory for distributed signal
//     shared_signal = (shared complex_t *) upc_all_alloc(THREADS, local_N * sizeof(complex_t));
//     if (shared_signal == NULL) {
//         fprintf(stderr, "Thread %d: Shared memory allocation error\n", MYTHREAD);
//         exit(1);
//     }

//     upc_barrier; // Wait for all allocations to complete

//     // Step 3: Distribute data from thread 0 to all threads
//     if (MYTHREAD == 0) {
//         printf("Thread 0: Distributing data to all threads\n");
//         for (int t = 0; t < THREADS; t++) {
//             for (int i = 0; i < local_N; i++) {
//                 shared_signal[t * local_N + i] = signal_full[t * local_N + i];
//             }
//         }
//         free(signal_full);
//         signal_full = NULL;
//         printf("Thread 0: Data distribution completed\n");
//     }

//     upc_barrier; // Wait for data distribution to complete

//     // Start timing
//     if (MYTHREAD == 0) {
//         gettimeofday(&tv_start, NULL);
//         printf("Thread 0: Starting FFT computation\n");
//     }

//     // Step 4: FFT computation
//     fft_compute_stages_upc(N, local_N);

//     // End timing
//     if (MYTHREAD == 0) {
//         gettimeofday(&tv_end, NULL);
//         start_time = tv_start.tv_sec + tv_start.tv_usec / 1000000.0;
//         end_time = tv_end.tv_sec + tv_end.tv_usec / 1000000.0;
//         printf("Thread 0: FFT computation completed\n");
//     }

//     upc_barrier;

//     // Step 5: Gather results back to thread 0
//     if (MYTHREAD == 0) {
//         printf("Thread 0: Gathering results\n");
//         for (int t = 0; t < THREADS; t++) {
//             for (int i = 0; i < local_N; i++) {
//                 fft_result[t * local_N + i] = shared_signal[t * local_N + i];
//             }
//         }
//     }

//     upc_barrier;

//     // Step 6: Display results (Thread 0)
//     if (MYTHREAD == 0) {
//         printf("Thread 0: FFT computation completed.\n");
//         printf("Thread 0: Execution time: %f seconds\n", end_time - start_time);

//         // Save result to file
//         FILE *outfile = fopen(OUTPUT, "w");
//         if (outfile) {
//             fprintf(outfile, "%d %d\n", N, Fs);
//             for (int i = 0; i < N; i++) {
//                 fprintf(outfile, "%.8f %.8f\n", fft_result[i].real, fft_result[i].imag);
//             }
//             fclose(outfile);
//             printf("Thread 0: Result saved to file '%s'.\n", OUTPUT);
//         } else {
//             perror("Cannot open output file");
//         }

//         free(fft_result);
//     }

//     // Step 7: Clean up shared memory
//     upc_free(shared_signal);

//     return 0;
// }
