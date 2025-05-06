import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_signal_and_fft(fft_data_file, signal_data_file, output_png):
    try:
        # Odczyt danych FFT
        with open(fft_data_file, 'r') as f_fft:
            header_line_fft = f_fft.readline().strip()
            n_fft, fs_fft = map(float, header_line_fft.split())
            n_fft = int(n_fft)
            fft_data = np.loadtxt(f_fft)
            real = fft_data[:, 0]
            imag = fft_data[:, 1]
            complex_spectrum = real + 1j * imag
            magnitude_spectrum = np.abs(complex_spectrum)
            phase_spectrum = np.arctan2(imag, real)

            frequencies = np.fft.fftfreq(n_fft, 1/fs_fft)
            positive_frequencies = frequencies[:n_fft//2]
            positive_magnitude = magnitude_spectrum[:n_fft//2]
            positive_phase = phase_spectrum[:n_fft//2]

        # Odczyt danych sygnału w czasie (pomijając pierwszą linię z Fs)
        try:
            with open(signal_data_file, 'r') as f_signal:
                fs_time_str = f_signal.readline().strip() # Odczytaj i zignoruj pierwszą linię
                time_signal = np.loadtxt(f_signal)
                if time_signal.ndim > 1:
                    time_signal = time_signal[:, 0] # Jeśli jest wiele kolumn, bierz pierwszą

                n_time = len(time_signal)
                time = np.linspace(0, n_time / fs_fft, n_time, endpoint=False) # Używamy Fs z FFT dla osi czasu
        except FileNotFoundError:
            print(f"Ostrzeżenie: Nie znaleziono pliku z sygnałem czasowym: {signal_data_file}")
            time = np.array([])
            time_signal = np.array([])
        except ValueError:
            print(f"Ostrzeżenie: Nieprawidłowy format pliku sygnału czasowego: {signal_data_file}.")
            time = np.array([])
            time_signal = np.array([])

        # Tworzenie podwykresów
        fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=False)

        # Wykres sygnału w dziedzinie czasu
        if time.size > 0 and time_signal.size > 0:
            axs[0].plot(time, time_signal)
            axs[0].set_xlabel("Czas (s)")
            axs[0].set_ylabel("Amplituda")
            axs[0].set_title("Sygnał w Dziedzinie Czasu")
            axs[0].grid(True)

        # Wykres widma amplitudowego
        axs[1].plot(positive_frequencies, positive_magnitude)
        axs[1].set_ylabel("Amplituda")
        axs[1].set_title("Widmo Amplitudowe FFT")
        axs[1].grid(True)

        # Wykres widma fazowego
        axs[2].plot(positive_frequencies, positive_phase)
        axs[2].set_xlabel("Częstotliwość (Hz)")
        axs[2].set_ylabel("Faza (radiany)")
        axs[2].grid(True)

        plt.tight_layout()
        plt.savefig(output_png)
        print(f"Wykres sygnału czasowego i widma zapisano do {output_png}")

    except FileNotFoundError:
        print(f"Błąd: Nie znaleziono pliku {fft_data_file}")
    except ValueError:
        print(f"Błąd: Nieprawidłowy format pliku {fft_data_file}. Pierwsza linia powinna zawierać 'N Fs'.")
    except Exception as e:
        print(f"Wystąpił błąd: {e}")

if __name__ == "__main__":
    if len(sys.argv) == 3:
        fft_data_file = sys.argv[1]
        signal_data_file = sys.argv[2]
        output_png = "fft_spectrum_with_time.png"
        plot_signal_and_fft(fft_data_file, signal_data_file, output_png)
    else:
        print("Użycie: python plot_fft.py <plik_z_danymi_fft> <plik_z_sygnalem_czasowym>")