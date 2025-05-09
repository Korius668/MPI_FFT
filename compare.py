import sys
import numpy as np
import matplotlib.pyplot as plt

def load_freq_amp_phase(filename):
    with open(filename) as f:
        header = f.readline()
        num_points, sampling_rate = map(float, header.strip().split())
        data = [tuple(map(float, line.split())) for line in f]
    return int(num_points), sampling_rate, data

def load_fft_data(filename):
    with open(filename) as f:
        header = f.readline()
        num_points, sampling_rate = map(float, header.strip().split())
        data = np.loadtxt(f)
    real, imag = data[:, 0], data[:, 1]
    return int(num_points), sampling_rate, real + 1j * imag

def main(file1, file2):
    from collections import defaultdict

    n1, sr1, data1 = load_freq_amp_phase(file1)
    freq_dict = defaultdict(list)
    for freq, amp, phase in data1:
        freq_dict[freq].append((amp, phase))

    freqs1 = np.fft.rfftfreq(n1, d=1.0/sr1)
    fft_bins = np.zeros_like(freqs1, dtype=complex)
    freq_to_index = {round(f, 8): i for i, f in enumerate(freqs1)}
    for freq, vals in freq_dict.items():
        freq_rounded = round(freq, 8)
        if freq_rounded in freq_to_index:
            idx = freq_to_index[freq_rounded]
            complex_sum = sum(amp * np.exp(1j * phase) for amp, phase in vals)
            fft_bins[idx] = complex_sum
    
    time_signal = np.fft.irfft(fft_bins, n=n1)
    time_axis = np.arange(n1) / sr1

    n2, sr2, fft_from_file = load_fft_data(file2)
    fft_from_file[n//2 + 1:] = fft_from_file[n//2 + 1:][::-1]
    reconstructed_signal = np.fft.ifft(fft_from_file).real
    t2 = np.linspace(0, n2 / sr2, n2, endpoint=False)
    
    plt.figure(figsize=(12, 8))
    
    plt.subplot(2, 2, 1)
    plt.plot(time_axis, time_signal)
    plt.title("Oryginalny sygnał")
    plt.xlabel("Czas [s]")
    plt.ylabel("Wartość")
    
    plt.subplot(2, 2, 2)
    magnitudes = np.abs(fft_bins)
    plt.bar(freqs1, magnitudes, width=freqs1[1] - freqs1[0])
    plt.title("Oryginalne Spektrum")
    plt.xlabel("Częstotliwość [Hz]")
    plt.ylabel("Amplituda")
    
    plt.subplot(2, 2, 3)
    plt.plot(t2, reconstructed_signal)
    plt.title("Zrekonstruowany sygnał")
    plt.xlabel("Czas [s]")
    plt.ylabel("Wartość")
    
    plt.subplot(2, 2, 4)
    amplitudes2 = np.sqrt(fft_from_file.real**2 + fft_from_file.imag**2)[:n2//2 + 1] / (n2 / 2)
    freqs2 = np.fft.rfftfreq(n2, d=1.0/sr2)
    plt.bar(freqs2, amplitudes2, width=freqs2[1] - freqs2[0])
    plt.title("Wyliczone Spektrum")
    plt.xlabel("Częstotliwość [Hz]")
    plt.ylabel("Amplituda")
    
    plt.tight_layout()
    plt.savefig("plot.png");
    


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) != 2:
        print("Użycie: python compare.py expected result")
        sys.exit(1)
    main(args[0], args[1])

