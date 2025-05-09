import sys
import numpy as np
import matplotlib.pyplot as plt

def load_freq_amp_phase(filename):
    with open(filename) as f:
        header = f.readline()
        num_points, sampling_rate = map(float, header.strip().split())
        data = np.loadtxt(f)
    freqs, amps, phases = data[:, 0], data[:, 1], data[:, 2]
    return int(num_points), sampling_rate, freqs, amps, phases

def load_fft_data(filename):
    with open(filename) as f:
        header = f.readline()
        num_points, sampling_rate = map(float, header.strip().split())
        data = np.loadtxt(f)
    real, imag = data[:, 0], data[:, 1]
    return int(num_points), sampling_rate, real + 1j * imag

def reconstruct_time_signal(freqs, amps, phases, num_points, sampling_rate):
    t = np.linspace(0, num_points / sampling_rate, num_points, endpoint=False)
    signal = np.zeros_like(t, dtype=np.complex128)
    for f, a, p in zip(freqs, amps, phases):
        # Construct complex exponential: A * exp(i*(2πft + φ))
        signal += 0.5 * a * np.exp(1j * (2 * np.pi * f * t + p))  # positive freq
        signal += 0.5 * a * np.exp(-1j * (2 * np.pi * f * t + p))  # negative freq (conj pair)
    return t, signal.real  # take real part


def main(file1, file2):

    n1, sr1, freqs, amps, phases = load_freq_amp_phase(file1)
    time, signal = reconstruct_time_signal(freqs, amps, phases, n1, sr1)

    fft_vals = np.fft.fft(signal)
    fft_freqs = np.fft.fftfreq(n1, d=1/sr1)
    
    n2, sr2, fft_from_file = load_fft_data(file2)
    reconstructed_signal = np.fft.ifft(fft_from_file).real
    t2 = np.linspace(0, n2 / sr2, n2, endpoint=False)
    
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(time, signal)
    plt.title("Oryginalny sygnał")
    plt.xlabel("Czas [s]")
    plt.ylabel("Wartość")
    
    plt.subplot(2, 2, 2)
    plt.bar(freqs, amps, width=freqs[1] - freqs[0])
    plt.title("Oryginalne Spektrum FFT")
    plt.xlabel("Częstotliwość [Hz]")
    plt.ylabel("Amplituda")
    
    plt.subplot(2, 2, 3)
    amplitudes2 = np.sqrt(fft_from_file.real**2 + fft_from_file.imag**2)[:n2//2 + 1]
    freqs2 = np.fft.rfftfreq(n2, d=1.0/sr2)
    plt.bar(freqs2, amplitudes2, width=freqs2[1] - freqs2[0])
    plt.title("Wyliczone Spektrum FFT")
    plt.xlabel("Częstotliwość [Hz]")
    plt.ylabel("Amplituda")
    
    plt.subplot(2, 2, 4)
    plt.plot(t2, reconstructed_signal)
    plt.title("Zrekonstruowany sygnał")
    plt.xlabel("Czas [s]")
    plt.ylabel("Wartość")
    
    plt.tight_layout()
    plt.savefig("plot.png");


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 2:
        print("Użycie: python compare.py expected result")
        sys.exit(1)
    main(args[0], args[1])

