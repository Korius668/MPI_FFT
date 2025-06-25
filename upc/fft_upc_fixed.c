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
shared complex_t *shared_signal;

int is_power_of_two(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Function to pad signal with zeros to nearest power of 2
complex_t* pad_to_power_of_two(complex_t *signal, int n, int *padded_n) {
    if (is_power_of_two(n)) {
        *padded_n = n;
        complex_t *padded_signal = (complex_t*)malloc(n * sizeof(complex_t));
        if (padded_signal == NULL) {
            fprintf(stderr, "Thread %d: Memory allocation error\n", MYTHREAD);
            exit(1);
        }
        for (int i = 0; i < n; i++) {
            padded_signal[i] = signal[i];
        }
        return padded_signal;
    } else {
        int next_power_of_two = 1;
        while (next_power_of_two < n) {
            next_power_of_two <<= 1;
        }
        *padded_n = next_power_of_two;
        complex_t *padded_signal = (complex_t*)malloc(next_power_of_two * sizeof(complex_t));
        if (padded_signal == NULL) {
            fprintf(stderr, "Thread %d: Memory allocation error\n", MYTHREAD);
            exit(1);
        }
        for (int i = 0; i < n; i++) {
            padded_signal[i] = signal[i];
        }
        for (int i = n; i < next_power_of_two; i++) {
            padded_signal[i] = make_complex(0.0, 0.0);
        }
        return padded_signal;
    }
}

// Helper function to swap two complex numbers
void swap(complex_t *a, complex_t *b) {
    complex_t temp = *a;
    *a = *b;
    *b = temp;
}

// Bit-reversal permutation (performed by thread 0)
void bit_reverse_permutation(complex_t *x, int N) {
    int bits = 0;
    while ((1 << bits) < N) {
        bits++;
    }

    for (int i = 0; i < N; i++) {
        int rev = 0;
        int temp = i;
        for (int j = 0; j < bits; j++) {
            rev = (rev << 1) | (temp & 1);
            temp >>= 1;
        }
        if (rev > i) {
            swap(&x[i], &x[rev]);
        }
    }
}

// POPRAWIONA funkcja FFT - kluczowe zmiany w komunikacji między wątkami
void fft_compute_stages_upc(int N, int local_N) {
    int stages = 0;
    while ((1 << stages) < N) {
        stages++;
    }

    complex_t W, Wm, t, u;
    int m, m2;

    // Iteracja przez etapy FFT
    for (int s = 1; s <= stages; s++) {
        m = 1 << s;  // Rozmiar aktualnego pod-FFT
        m2 = m >> 1; // Połowa rozmiaru
        Wm = complex_exp(-2.0 * M_PI / m); // Współczynnik obrotowy dla etapu

        // Iteracja przez motylki w etapie
        for (int k = 0; k < N; k += m) {
            W = make_complex(1.0, 0.0); // Reset współczynnika obrotowego
            for (int j = 0; j < m2; j++) {
                int idx1 = k + j;
                int idx2 = k + j + m2;

                // Określenie, które wątki posiadają potrzebne indeksy
                int thread1 = idx1 / local_N;
                int thread2 = idx2 / local_N;
                int local_idx1 = idx1 % local_N;
                int local_idx2 = idx2 % local_N;

                if (thread1 == thread2) {
                    // Oba indeksy są na tym samym wątku
                    if (MYTHREAD == thread1) {
                        t = complex_mul(W, shared_signal[MYTHREAD * local_N + local_idx2]);
                        u = shared_signal[MYTHREAD * local_N + local_idx1];
                        shared_signal[MYTHREAD * local_N + local_idx1] = complex_add(u, t);
                        shared_signal[MYTHREAD * local_N + local_idx2] = complex_sub(u, t);
                    }
                } else {
                    // Indeksy są na różnych wątkach - POPRAWIONA LOGIKA
                    if (MYTHREAD == thread1) {
                        // Ten wątek ma idx1, potrzebuje idx2 z thread2
                        complex_t val_idx1 = shared_signal[MYTHREAD * local_N + local_idx1];
                        complex_t val_idx2 = shared_signal[thread2 * local_N + local_idx2];
                        
                        // Oblicz nową wartość dla idx1
                        t = complex_mul(W, val_idx2);
                        u = val_idx1;
                        shared_signal[MYTHREAD * local_N + local_idx1] = complex_add(u, t);
                    }
                    
                    if (MYTHREAD == thread2) {
                        // Ten wątek ma idx2, potrzebuje idx1 z thread1
                        complex_t val_idx1 = shared_signal[thread1 * local_N + local_idx1];
                        complex_t val_idx2 = shared_signal[MYTHREAD * local_N + local_idx2];
                        
                        // Oblicz nową wartość dla idx2
                        t = complex_mul(W, val_idx2);
                        u = val_idx1;
                        shared_signal[MYTHREAD * local_N + local_idx2] = complex_sub(u, t);
                    }
                }
                
                // Aktualizuj współczynnik obrotowy
                W = complex_mul(W, Wm);
            }
        }
        
        // KRYTYCZNA synchronizacja po każdym etapie
        upc_barrier;
    }
}

int main(int argc, char *argv[]) {
    int local_N;
    complex_t *signal_full = NULL;
    complex_t *fft_result = NULL;
    char *filename = NULL;
    FILE *infile = NULL;
    struct timeval tv_start, tv_end;
    double start_time, end_time;
    int temp_N = 0;
    
    if (MYTHREAD == 0) {
        printf("Thread 0: Starting UPC FFT with %d threads\n", THREADS);
    }
    
    // Krok 1: Odczyt danych (Wątek 0)
    if (MYTHREAD == 0) {
        if (argc != 2) {
            fprintf(stderr, "Usage: %s <signal_file>\n", argv[0]);
            exit(1);
        }
        filename = argv[1];
        infile = fopen(filename, "r");
        if (infile == NULL) {
            perror("Cannot open input file");
            exit(1);
        }

        int temp_Fs;
        if (fscanf(infile, "%d", &temp_Fs) != 1) {
            fprintf(stderr, "Error reading sampling frequency from file.\n");
            fclose(infile);
            exit(1);
        }
        Fs = temp_Fs;
        printf("Thread 0: Read sampling frequency: %d Hz.\n", Fs);

        // Odczytaj dane
        double real_val;
        int capacity = 1024;
        signal_full = (complex_t *)malloc(capacity * sizeof(complex_t));
        if (!signal_full) {
            perror("Memory allocation error (signal_full)");
            fclose(infile);
            exit(1);
        }

        while (fscanf(infile, "%lf", &real_val) == 1) {
            if (temp_N >= capacity) {
                capacity *= 2;
                complex_t *temp = (complex_t *)realloc(signal_full, capacity * sizeof(complex_t));
                if (!temp) {
                    perror("Memory reallocation error (signal_full)");
                    free(signal_full);
                    fclose(infile);
                    exit(1);
                }
                signal_full = temp;
            }
            signal_full[temp_N] = make_complex(real_val, 0.0);
            temp_N++;
        }
        fclose(infile);

        if (temp_N == 0) {
            fprintf(stderr, "Error: Input file is empty.\n");
            free(signal_full);
            exit(1);
        }
        
        // Dopełnij do potęgi 2
        int padded_N;
        complex_t *padded_signal_full = pad_to_power_of_two(signal_full, temp_N, &padded_N);
        if (padded_signal_full != signal_full) {
            free(signal_full);
            signal_full = padded_signal_full;
            temp_N = padded_N;
            printf("Thread 0: Signal padded with zeros to length %d\n", temp_N);
        }

        // Sprawdź podzielność przez liczbę wątków
        if (temp_N % THREADS != 0) {
            fprintf(stderr, "Error: Number of samples (%d) must be divisible by number of threads (%d).\n", temp_N, THREADS);
            free(signal_full);
            exit(1);
        }

        N = temp_N;
        printf("Thread 0: Processing %d samples with %d threads\n", N, THREADS);

        // WAŻNE: Bit-reversal permutation PRZED dystrybucją
        bit_reverse_permutation(signal_full, N);
        printf("Thread 0: Bit-reversal permutation completed\n");

        fft_result = (complex_t *)malloc(N * sizeof(complex_t));
        if (!fft_result) {
            perror("Memory allocation error (fft_result)");
            free(signal_full);
            exit(1);
        }
    }

    upc_barrier;

    local_N = N / THREADS;
    
    if (MYTHREAD == 0) {
        printf("Thread 0: Local N per thread: %d\n", local_N);
    }

    // Alokuj pamięć współdzieloną
    shared_signal = (shared complex_t *) upc_all_alloc(THREADS, local_N * sizeof(complex_t));
    if (shared_signal == NULL) {
        fprintf(stderr, "Thread %d: Shared memory allocation error\n", MYTHREAD);
        exit(1);
    }

    upc_barrier;

    // Dystrybuuj dane
    if (MYTHREAD == 0) {
        printf("Thread 0: Distributing data\n");
        for (int t = 0; t < THREADS; t++) {
            for (int i = 0; i < local_N; i++) {
                shared_signal[t * local_N + i] = signal_full[t * local_N + i];
            }
        }
        free(signal_full);
        signal_full = NULL;
    }

    upc_barrier;

    // Pomiar czasu
    if (MYTHREAD == 0) {
        gettimeofday(&tv_start, NULL);
        printf("Thread 0: Starting FFT computation\n");
    }

    // GŁÓWNE OBLICZENIA FFT
    fft_compute_stages_upc(N, local_N);

    if (MYTHREAD == 0) {
        gettimeofday(&tv_end, NULL);
        start_time = tv_start.tv_sec + tv_start.tv_usec / 1000000.0;
        end_time = tv_end.tv_sec + tv_end.tv_usec / 1000000.0;
        printf("Thread 0: FFT computation completed\n");
    }

    upc_barrier;

    // Zbierz wyniki
    if (MYTHREAD == 0) {
        for (int t = 0; t < THREADS; t++) {
            for (int i = 0; i < local_N; i++) {
                fft_result[t * local_N + i] = shared_signal[t * local_N + i];
            }
        }
    }

    upc_barrier;

    // Wyświetl wyniki
    if (MYTHREAD == 0) {
        printf("Thread 0: FFT computation completed.\n");
        printf("Thread 0: Execution time: %f seconds\n", end_time - start_time);

        FILE *outfile = fopen(OUTPUT, "w");
        if (outfile) {
            fprintf(outfile, "%d %d\n", N, Fs);
            for (int i = 0; i < N; i++) {
                fprintf(outfile, "%.8f %.8f\n", fft_result[i].real, fft_result[i].imag);
            }
            fclose(outfile);
            printf("Thread 0: Result saved to file '%s'.\n", OUTPUT);
        } else {
            perror("Cannot open output file");
        }

        free(fft_result);
    }

    upc_free(shared_signal);
    return 0;
}
