#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h> // Dla obsługi liczb zespolonych (double complex)
#include <mpi.h>     // Dla MPICH

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef OUTPUT
#define	OUTPUT "fft_output"
#endif

int is_power_of_two(int n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Funkcja dopełniająca sygnał zerami do najbliższej potęgi 2
double complex* pad_to_power_of_two(double complex *signal, int n, int *padded_n) {
    if (is_power_of_two(n)) {
        *padded_n = n;
        double complex *padded_signal = (double complex*)malloc(n * sizeof(double complex));
        if (padded_signal == NULL) {
            perror("Błąd alokacji pamięci");
            MPI_Abort(MPI_COMM_WORLD, 1);
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
        double complex *padded_signal = (double complex*)malloc(next_power_of_two * sizeof(double complex));
        if (padded_signal == NULL) {
            perror("Błąd alokacji pamięci");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for (int i = 0; i < n; i++) {
            padded_signal[i] = signal[i];
        }
        for (int i = n; i < next_power_of_two; i++) {
            padded_signal[i] = 0.0 + 0.0 * I;
        }
        return padded_signal;
    }
}

// Funkcja pomocnicza do zamiany dwóch liczb zespolonych
void swap(double complex *a, double complex *b) {
    double complex temp = *a;
    *a = *b;
    *b = temp;
}

// Funkcja wykonująca permutację bit-reversal na całym sygnale (wykonywane przez proces 0)
void bit_reverse_permutation(double complex *x, int N) {
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

// Funkcja obliczająca część iteracyjnego algorytmu FFT (motylki)
// Działa na lokalnych danych, ale wymaga komunikacji dla motylków międzyprocesowych
void fft_compute_stages(double complex *local_x, int N, int local_N, int rank, int size) {
    int stages = 0;
    while ((1 << stages) < N) {
        stages++;
    }

    double complex W, Wm, t, u;
    int m, m2;
    MPI_Status status;

    // Iteracja przez etapy FFT
    for (int s = 1; s <= stages; s++) {
        m = 1 << s;  // Rozmiar aktualnego pod-FFT
        m2 = m >> 1; // Połowa rozmiaru
        Wm = cexp(-2.0 * M_PI * I / m); // Współczynnik obrotowy dla etapu

        // Iteracja przez motylki w etapie
        for (int k = 0; k < N; k += m) {
            W = 1.0 + 0.0 * I; // Reset współczynnika obrotowego dla każdej grupy motylków
            for (int j = 0; j < m2; j++) {
                int idx1 = k + j;
                int idx2 = k + j + m2;

                // Określenie, które procesy posiadają potrzebne indeksy
                int rank1 = idx1 / local_N;
                int rank2 = idx2 / local_N;
                int local_idx1 = idx1 % local_N;
                int local_idx2 = idx2 % local_N;

                if (rank == rank1 && rank == rank2) {
                    // Oba indeksy są lokalne - standardowy motylek
                    t = W * local_x[local_idx2];
                    u = local_x[local_idx1];
                    local_x[local_idx1] = u + t;
                    local_x[local_idx2] = u - t;
                } else if (rank == rank1) {
                    // Ten proces (rank) ma idx1, potrzebuje idx2 od rank2
                    double complex val_idx2;
                    // Wyślij moją wartość idx1 do partnera (rank2)
                    MPI_Sendrecv(&local_x[local_idx1], 1, MPI_DOUBLE_COMPLEX, rank2, 0,
                                 &val_idx2,          1, MPI_DOUBLE_COMPLEX, rank2, 0,
                                 MPI_COMM_WORLD, &status);

                    // Oblicz moją część motylka używając otrzymanej wartości
                    t = W * val_idx2;
                    u = local_x[local_idx1]; // Moja oryginalna wartość
                    local_x[local_idx1] = u + t; // Aktualizuj tylko moją część
                } else if (rank == rank2) {
                     // Ten proces (rank) ma idx2, potrzebuje idx1 od rank1
                    double complex val_idx1;
                    // Odbierz wartość idx1 od partnera (rank1) i wyślij mu moją idx2
                    // Uwaga: Kolejność Send/Recv musi być spójna z partnerem!
                    MPI_Sendrecv(&local_x[local_idx2], 1, MPI_DOUBLE_COMPLEX, rank1, 0,
                                 &val_idx1,          1, MPI_DOUBLE_COMPLEX, rank1, 0,
                                 MPI_COMM_WORLD, &status);

                    // Oblicz moją część motylka używając otrzymanej wartości
                    t = W * local_x[local_idx2]; // Moja oryginalna wartość
                    u = val_idx1;                // Wartość otrzymana od partnera
                    local_x[local_idx2] = u - t; // Aktualizuj tylko moją część
                }
                 // Aktualizuj współczynnik obrotowy dla następnego motylka w grupie
                 W *= Wm;
            }
        }
         // Bariera synchronizacyjna po każdym etapie jest ważna,
         // aby upewnić się, że wszystkie komunikacje Sendrecv zostały zakończone
         // przed przejściem do następnego etapu.
         MPI_Barrier(MPI_COMM_WORLD);
    }
}


int main(int argc, char *argv[]) {
    int rank, size;
    int N = 0; // Całkowity rozmiar sygnału
    int local_N; // Rozmiar sygnału na lokalnym procesie
    double complex *signal_full = NULL; // Pełny sygnał (tylko w procesie 0)
    double complex *local_signal = NULL; // Lokalna część sygnału
    double complex *fft_result = NULL;   // Wynik FFT (tylko w procesie 0)
    char *filename = NULL;
    FILE *infile = NULL;
    double start_time, end_time;
    int Fs;
    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // --- Krok 1: Odczyt danych i dystrybucja (Proces 0) ---
    if (rank == 0) {
        if (argc != 2) {
            fprintf(stderr, "Użycie: mpiexec -np <liczba_procesów> %s <plik_z_sygnalem>\n", argv[0]);
            MPI_Abort(MPI_COMM_WORLD, 1); // Zakończ wszystkie procesy
        }
        filename = argv[1];
        infile = fopen(filename, "r");
        if (infile == NULL) {
            perror("Nie można otworzyć pliku wejściowego");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        if (fscanf(infile, "%d", &Fs) != 1) {
            fprintf(stderr, "Błąd odczytu częstotliwości próbkowania z pliku.\n");
            fclose(infile);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Proces 0: Odczytano częstotliwość próbkowania: %d Hz.\n", Fs);

        // Odczytaj dane i policz N
        double real_val;
        int capacity = 1024;
        signal_full = (double complex *)malloc(capacity * sizeof(double complex));
        if (!signal_full) {
             perror("Błąd alokacji pamięci (signal_full)");
             fclose(infile);
             MPI_Abort(MPI_COMM_WORLD, 1);
        }

        while (fscanf(infile, "%lf", &real_val) == 1) {
            if (N >= capacity) {
                capacity *= 2;
                double complex *temp = (double complex *)realloc(signal_full, capacity * sizeof(double complex));
                 if (!temp) {
                    perror("Błąd realokacji pamięci (signal_full)");
                    free(signal_full);
                    fclose(infile);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                 }
                 signal_full = temp;
            }
            signal_full[N] = real_val + 0.0 * I; // Zapisz jako liczbę zespoloną
            N++;
        }
        fclose(infile);

        // Sprawdź warunki N
        if (N == 0) {
             fprintf(stderr, "Błąd: Plik wejściowy jest pusty.\n");
             free(signal_full);
             MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Sprawdzenie czy N jest potęgą 2
        int padded_N;
        double complex *padded_signal_full = pad_to_power_of_two(signal_full, N, &padded_N);
        if (padded_signal_full != signal_full) { // Jeśli nastąpiło dopełnienie
            free(signal_full);
            signal_full = padded_signal_full;
            N = padded_N;
            printf("Proces 0: Sygnał dopełniono zerami do długości %d (najbliższa potęga 2).\n", N);
        } else {
            printf("Proces 0: Liczba próbek (%d) jest już potęgą 2.\n", N);
        }

        // Sprawdzenie czy N jest podzielne przez liczbę procesów
        if (N % size != 0) {
            fprintf(stderr, "Błąd: Liczba próbek (%d) musi być podzielna przez liczbę procesów (%d).\n", N, size);
            free(signal_full);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        printf("Proces 0: Odczytano %d próbek z pliku '%s'. Liczba procesów: %d.\n", N, filename, size);

        // Wykonaj permutację bit-reversal na pełnym sygnale przed rozproszeniem
        bit_reverse_permutation(signal_full, N);

        // Alokuj pamięć na wynik końcowy
        fft_result = (double complex *)malloc(N * sizeof(double complex));
         if (!fft_result) {
             perror("Błąd alokacji pamięci (fft_result)");
             free(signal_full);
             MPI_Abort(MPI_COMM_WORLD, 1);
         }

    } // koniec if (rank == 0)

    start_time = MPI_Wtime(); // Rozpocznij pomiar czasu

    // --- Krok 2: Rozgłoś rozmiar N do wszystkich procesów ---
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Oblicz lokalny rozmiar
    local_N = N / size;

    // --- Krok 3: Alokuj lokalne bufory i rozprosz dane ---
    local_signal = (double complex *)malloc(local_N * sizeof(double complex));
    if (!local_signal) {
        perror("Błąd alokacji pamięci (local_signal)");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Rozprosz dane (po permutacji bit-reversal) z procesu 0 do wszystkich procesów
    MPI_Scatter(signal_full, local_N, MPI_DOUBLE_COMPLEX, // Dane wysyłane
                local_signal, local_N, MPI_DOUBLE_COMPLEX, // Bufor odbiorczy
                0, MPI_COMM_WORLD);

    // Proces 0 może teraz zwolnić pamięć pełnego sygnału wejściowego
    if (rank == 0) {
        free(signal_full);
        signal_full = NULL; // Dobra praktyka
    }

    // --- Krok 4: Obliczenia FFT ---
    // Każdy proces wykonuje obliczenia na swoich danych, komunikując się w razie potrzeby
    fft_compute_stages(local_signal, N, local_N, rank, size);


    // --- Krok 5: Zbierz wyniki z powrotem do procesu 0 ---
    MPI_Gather(local_signal, local_N, MPI_DOUBLE_COMPLEX, // Dane lokalne do wysłania
               fft_result,   local_N, MPI_DOUBLE_COMPLEX, // Bufor zbiorczy (tylko rank 0)
               0, MPI_COMM_WORLD);


    end_time = MPI_Wtime(); // Zakończ pomiar czasu

    // --- Krok 6: Wyświetl wyniki (Proces 0) ---
    if (rank == 0) {
        printf("Proces 0: Zakończono obliczenia FFT.\n");
        printf("Proces 0: Czas wykonania: %f sekund\n", end_time - start_time);

        // Opcjonalnie: Wypisz wynik FFT
        // printf("Wynik FFT:\n");
        // for (int i = 0; i < N; i++) {
        //     printf("[%d]: %.6f + %.6fi\n", i, creal(fft_result[i]), cimag(fft_result[i]));
        // }

        // Opcjonalnie: Zapisz wynik do pliku
        FILE *outfile = fopen(OUTPUT, "w");
        if (outfile) {
             fprintf(outfile, "%d %d\n", N, Fs);
             for (int i = 0; i < N; i++) {
                 fprintf(outfile, "%.8f %.8f\n", creal(fft_result[i]), cimag(fft_result[i]));
             }
             fclose(outfile);
             printf("Proces 0: Wynik zapisano do pliku 'fft_output'.\n");
        } else {
            perror("Nie można otworzyć pliku wyjściowego 'fft_output'");
        }

        free(fft_result); // Zwolnij pamięć wyniku
    }

    // --- Krok 7: Zwolnij pamięć lokalną i zakończ MPI ---
    free(local_signal);
    MPI_Finalize();

    return 0;
}
