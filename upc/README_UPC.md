# UPC_FFT
Program do obliczania FFT w środowisku rozproszonym z użyciem biblioteki UPC (Unified Parallel C) dla programowania PGAS w laboratorium 204.

## O UPC
UPC (Unified Parallel C) to język programowania równoległego wykorzystujący model PGAS (Partitioned Global Address Space). W przeciwieństwie do podejścia MPI opartego na przekazywaniu komunikatów, UPC pozwala wątkom na bezpośredni dostęp do pamięci na zdalnych wątkach poprzez współdzieloną przestrzeń adresową.

## Kluczowe różnice w stosunku do wersji MPI
- **Globalna przestrzeń adresowa**: Bezpośredni dostęp do pamięci zamiast przekazywania komunikatów
- **Tablice współdzielone**: Dane rozproszone między wątkami używając słowa kluczowego `shared`
- **Uproszczona komunikacja**: Brak jawnych operacji send/receive
- **Wbudowana synchronizacja**: `upc_barrier` do synchronizacji wątków
- **Identyfikacja wątków**: `MYTHREAD` i `THREADS` zamiast rank i size

## Przygotowanie środowiska w Lab 204

### 1. Załadowanie zmiennych środowiskowych
\`\`\`bash
# Opcjonalnie, dla współpracy z MPI
source /opt/nfs/config/source_mpich411.sh

# Główne środowisko Berkeley UPC
source /opt/nfs/config/source_bupc.sh
\`\`\`

### 2. Przygotowanie pliku z listą węzłów
\`\`\`bash
# Automatyczne generowanie (przez Makefile)
make nodes

# Lub ręcznie
/opt/nfs/config/station204_name_list.sh 1 16 > nodes
\`\`\`

### 3. Sprawdzenie środowiska
\`\`\`bash
make check-env
\`\`\`

## Generowanie przykładowych danych wejściowych
\`\`\`bash
make generate
\`\`\`

To utworzy pliki example1, example2 i example3 z różnymi charakterystykami sygnałów.

## Obsługa Makefile

### Kompilacja
\`\`\`bash
# Kompilacja z Berkeley UPC (domyślnie)
make compile-bupc

# Kompilacja z GNU UPC (alternatywnie)
make compile-gupc
\`\`\`

### Uruchomienie
\`\`\`bash
# Uruchomienie z domyślnymi parametrami (4 węzły × 4 wątki = 16 wątków)
make run

# Niestandardowa liczba węzłów
make run NODES=8

# Niestandardowa liczba wątków na proces
make run PTHREADS=2

# Niestandardowy plik wejściowy
make run INPUT=example2

# Kombinacja parametrów
make run NODES=6 PTHREADS=4 INPUT=example3
\`\`\`

### Testowanie różnych konfiguracji
\`\`\`bash
# Mała konfiguracja: 4 wątki (2 węzły × 2 wątki/proces)
make test-small

# Średnia konfiguracja: 16 wątków (4 węzły × 4 wątki/proces)
make test-medium

# Duża konfiguracja: 32 wątków (8 węzłów × 4 wątki/proces)
make test-large
\`\`\`

### Porównanie wyników
\`\`\`bash
make compare
\`\`\`

### Czyszczenie
\`\`\`bash
# Usunięcie skompilowanych plików
make clean

# Przywrócenie stanu początkowego
make reset
\`\`\`

## Szczegóły kompilacji i uruchomienia

### Kompilacja Berkeley UPC
\`\`\`bash
upcc -bupc -network=udp -pthreads=4 fft_upc.c -o fft_bupc -lm
\`\`\`

### Kompilacja GNU UPC
\`\`\`bash
/opt/nfs/berkeley_upc-2021.4.0/bin/upcc -gupc -Wc,"-fPIE" -network=udp -pthreads=4 fft_upc.c -o fft_gupc -lm
\`\`\`

### Uruchomienie
\`\`\`bash
UPC_NODEFILE=nodes upcrun -shared-heap 256M -c 4 -N 4 -n 16 ./fft_bupc example1
\`\`\`

Gdzie:
- `-shared-heap 256M`: Rozmiar współdzielonej sterty
- `-c 4`: Liczba wątków UPC na proces
- `-N 4`: Liczba węzłów
- `-n 16`: Całkowita liczba wątków (4×4=16)

## Model programowania UPC

### Model pamięci
- **Pamięć prywatna**: Każdy wątek ma własną prywatną przestrzeń pamięci
- **Pamięć współdzielona**: Rozproszona między wszystkimi wątkami, dostępna dla każdego wątku
- **Koligacja**: Każdy współdzielony obiekt ma koligację z określonym wątkiem

### Kluczowe funkcje UPC użyte w programie
1. **Tablice współdzielone**: `shared [] double complex *shared_signal`
2. **Identyfikacja wątków**: `MYTHREAD` (bieżący wątek), `THREADS` (całkowita liczba wątków)
3. **Alokacja pamięci**: `upc_all_alloc()` dla pamięci współdzielonej
4. **Synchronizacja**: `upc_barrier` dla synchronizacji wątków
5. **Pomiar czasu**: `upc_wtime()` dla pomiaru wydajności

### Wzorzec komunikacji
Wersja UPC eliminuje jawne przekazywanie komunikatów poprzez:
- Użycie tablic współdzielonych do dystrybucji danych
- Bezpośredni dostęp do pamięci dla operacji motylkowych
- Synchronizację barierową zamiast komunikacji punkt-punkt

## Rozważania dotyczące wydajności
- **Lokalność pamięci**: UPC automatycznie obsługuje umieszczanie danych
- **Narzut komunikacyjny**: Zmniejszony w porównaniu z MPI dzięki bezpośredniemu dostępowi do pamięci
- **Skalowalność**: Lepiej dostosowany do systemów z pamięcią współdzieloną lub architektur NUMA

## Wymagania kompilacji
- Kompilator UPC (Berkeley UPC lub GNU UPC)
- Wsparcie biblioteki matematycznej dla liczb zespolonych
- Wsparcie standardu C99

## Przykłady użycia

### Podstawowe użycie
\`\`\`bash
# Przygotowanie środowiska
source /opt/nfs/config/source_bupc.sh
make setup
make nodes

# Generowanie danych testowych
make generate

# Kompilacja i uruchomienie
make compile-bupc
make run

# Porównanie wyników
make compare
\`\`\`

### Testowanie wydajności
\`\`\`bash
# Test z różnymi konfiguracjami
make test-small    # 4 wątki
make test-medium   # 16 wątków  
make test-large    # 32 wątki

# Niestandardowa konfiguracja
make run NODES=12 PTHREADS=4 INPUT=example2  # 48 wątków
\`\`\`

### Rozwiązywanie problemów
\`\`\`bash
# Sprawdzenie środowiska
make check-env

# Sprawdzenie dostępnych węzłów
cat nodes

# Uruchomienie z mniejszą stertą współdzieloną
make run SHARED_HEAP=128M
\`\`\`

## Struktura plików
- `fft_upc.c` - Główny program UPC FFT
- `sygnalGenerate.py` - Generator przykładowych sygnałów
- `compare.py` - Skrypt porównujący wyniki
- `Makefile` - Automatyzacja kompilacji i uruchomienia
- `nodes` - Lista dostępnych węzłów (generowana automatycznie)

## Wsparcie
W przypadku problemów sprawdź:
1. Czy środowisko UPC jest prawidłowo załadowane (`make check-env`)
2. Czy plik węzłów istnieje (`make nodes`)
3. Czy parametry są prawidłowe (liczba wątków musi być podzielna przez liczbę próbek)
4. Czy rozmiar współdzielonej sterty jest wystarczający
