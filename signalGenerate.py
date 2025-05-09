import numpy as np
import sys

def generate_sine_wave(freq, sample_rate, duration, amplitude=1.0, phase=0):
    t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    return amplitude * np.sin(2 * np.pi * freq * t + phase)

def save_signal_to_file(filename, sample_rate, signal):
    with open(filename, 'w') as f:
        f.write(f"{sample_rate}\n")  
        for sample in signal:
            f.write(f"{sample}\n")

def generate(n, filename):
    m = n // 2
    freqs = np.random.randint(1, m, size=n)
    phases = np.random.rand(n) *  np.pi
    amplitudes = np.random.rand(n) * 10
    wave = sum(generate_sine_wave(freq, n, 1.0, amp, phase) for freq, amp, phase in zip(freqs, amplitudes, phases))
    fft_res = np.fft.fft(wave)
    reals = fft_res.real
    imags = fft_res.imag
    with open(filename + "_expected", "w") as file:
        file.write(f"{n} {n}\n")
        file.write('\n'.join(map(lambda x: ' '.join(map(str, x)), zip(freqs, amplitudes, phases, reals, imags)) ))
    save_signal_to_file(filename, n, wave)
    

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) == 0:
        main()
    elif len(args) == 2:
        try:
            n = int(args[0])
            assert n >= 4
        except:
            print("Error: n musi być liczbą całkowitą >= 4")
            sys.exit(1)

        filename = 'generated_' + args[1]
        generate(n, filename)
    else:
        print("Error: niewłaściwa liczba argumentów")
        sys.exit(1)
