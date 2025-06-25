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
// Main FFT array and a dummy one used for alternating read/write operations
shared [BLOCKSIZE] complex_t *shared_signal;
shared [BLOCKSIZE] complex_t *tmp_ptr;

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

void fft_compute_stages_upc(int N, shared [BLOCKSIZE] complex_t *x, shared [BLOCKSIZE] complex_t *tmp) {
    int stages = 0;
    while ((1 << stages) < N) stages++;

    shared[BLOCKSIZE] complex_t *read_buf = x;
    shared[BLOCKSIZE] complex_t *write_buf = tmp;

    for (int s = 1; s <= stages; s++) {
        int m = 1 << s;
        int m2 = m >> 1;
        complex_t Wm = complex_exp(-2.0 * M_PI / m);

        // Nice substitution for checking (MYTHREAD == ...)
	upc_forall(int k = 0; k < N; k += m; &read_buf[k]) {
            complex_t W = make_complex(1.0, 0.0);
            for (int j = 0; j < m2; j++) {
                int idx1 = k + j;
                int idx2 = k + j + m2;

		if (idx2 >= N) continue;

                complex_t u = read_buf[idx1];
                complex_t v = read_buf[idx2];
                complex_t t = complex_mul(W, v);

                write_buf[idx1] = complex_add(u, t);
                write_buf[idx2] = complex_sub(u, t);

            	W = complex_mul(W, Wm);
            }
        }
        upc_barrier;

        // Alternating read/write buffers to avoid race conditions
	shared [BLOCKSIZE] complex_t *p_tmp = read_buf;
	read_buf = write_buf;
	write_buf = p_tmp;
    }
    if (stages % 2 == 1) {
        upc_forall(int i = 0; i < N; i++; &x[i]) {
		x[i] = read_buf[i];
	}
	upc_barrier;
    }
}

int main(int argc, char *argv[]) {
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
        fft_result = malloc(N * sizeof(complex_t));
    }

    upc_barrier;

    tmp_ptr = (shared [BLOCKSIZE] complex_t *) upc_all_alloc(1, N * sizeof(complex_t));
    shared_signal = (shared [BLOCKSIZE] complex_t *) upc_all_alloc(1, N * sizeof(complex_t));

    if (MYTHREAD == 0) {
	    if (tmp_ptr == NULL || shared_signal == NULL) {
	        fprintf(stderr, "Thread %d: shared_signal allocation failed\n", MYTHREAD);
	        upc_global_exit(1);
	    }
    }

    upc_barrier;

    // Bit-Reverse Permutation
    if (MYTHREAD == 0) {
	int bits = 0;
	while ((1 << bits) < N) bits++;
	for (int i = 0; i < N; i++) {
		int rev = 0, tmp = i;
		for (int j = 0; j < bits; j++) {
			rev = (rev << 1) | (tmp & 1);
			tmp >>= 1;
		}
		shared_signal[rev] = signal_full[i];
	}
        free(signal_full);
    }
    upc_barrier;

    // FFT Execution
    if (MYTHREAD == 0) gettimeofday(&tv_start, NULL);
    fft_compute_stages_upc(N, shared_signal, tmp_ptr);
    if (MYTHREAD == 0) gettimeofday(&tv_end, NULL);

    upc_barrier;
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

    upc_all_free(shared_signal);
    upc_all_free(tmp_ptr);
    return 0;
}
