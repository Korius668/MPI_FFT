# MPI_FFT
Program do obliczania FFT w środowisku na środowisku rozproszonym z użyciem biblioteki OpenMPI.
## Połączenie się z maszyną w sali i przygotowanie środowiska

Aby połączyć się z maszyną w sali, użyj polecenia `ssh`:
```bash
ssh <nazwa_komputera_w_sali>
```
Na przykład:
```bash
ssh stud204-09
```
Następnie, aby przygotować środowisko, wykonaj następujące polecenia:

```bash
source /opt/nfs/config/source_mpich430.sh
source /opt/nfs/config/source_cuda121.sh
```
# Generowanie przykładowych danych wejściowych
Aby wygenerować przykładowe dane wejściowe, użyj programu sygnalGenerate.py.

```bash
python sygnalGenerate.py
```
Jeśli zależy Ci na innych sygnałach, zmodyfikuj kod tego programu.

# Obsługa Makefile
Kompilacja
Aby skompilować projekt, użyj polecenia:
```bash
make compile
```
Uruchomienie
Aby uruchomić program z domyślnym inputem, użyj polecenia:
```bash
make run
```
Dla niestandardowego inputu, użyj:
```bash
make run INPUT=<nazwa_pliku>
```
Na przykład:
```bash
make run INPUT=example_signals/sinus
```
Czyszczenie
Aby wyczyścić miejsce robocze (usunąć skompilowane pliki), wykonaj:
```bash
make clean
```
Resetowanie
Aby przywrócić miejsce robocze do stanu początkowego (usunąć również wygenerowane dane), wykonaj:
```bash
make reset
```
