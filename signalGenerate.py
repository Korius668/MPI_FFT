import numpy as np

def generate_sine_wave(freq, sample_rate, duration, amplitude=1.0, phase=0):
    t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    return amplitude * np.sin(2 * np.pi * freq * t + phase)

def save_signal_to_file(filename, sample_rate, signal):
    with open(filename, 'w') as f:
        f.write(f"{sample_rate}\n")  
        for sample in signal:
            f.write(f"{sample}\n")

def main():
    sample_rate = 44100  # Hz
    duration = 1.0       # sekundy


    freq1 = 440          # Hz
    freq2 = 880          # Hz


    sine1 = generate_sine_wave(freq1, sample_rate, duration)
    sine2 = generate_sine_wave(freq2, sample_rate, duration)
    sum_sines = sine1 + sine2
    modulated = sine1 * sine2
    composite = sine1 + 0.5 * sine2 + 0.25 * generate_sine_wave(1760, sample_rate, duration)

    save_signal_to_file("sine1.txt", sample_rate, sine1)
    save_signal_to_file("sine2.txt", sample_rate, sine2)
    save_signal_to_file("sum_sines.txt", sample_rate, sum_sines)
    save_signal_to_file("modulated.txt", sample_rate, modulated)
    save_signal_to_file("composite.txt", sample_rate, composite)

if __name__ == "__main__":
    main()
