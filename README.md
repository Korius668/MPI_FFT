# Połączenie się z maszyną w sali i przygotowanie środowiska

Aby połączyć się z maszyną w sali, użyj polecenia `ssh`:

```bash
ssh <nazwa_komputera_w_sali>
Na przykład:

Bash

ssh stud204-09
Następnie, aby przygotować środowisko, wykonaj następujące polecenia:

Bash

source /opt/nfs/config/source_mpich430.sh
source /opt/nfs/config/source_cuda121.sh
Generowanie przykładowych danych wejściowych
Aby wygenerować przykładowe dane wejściowe, użyj programu sygnalGenerate.py.

Bash

python sygnalGenerate.py
Jeśli zależy Ci na innych sygnałach, zmodyfikuj kod tego programu.

Obsługa Makefile
Kompilacja
Aby skompilować projekt, użyj polecenia:

Bash

make compile
Uruchomienie
Aby uruchomić program z domyślnym inputem, użyj polecenia:

Bash

make run
Dla niestandardowego inputu, użyj:

Bash

make run INPUT=<nazwa_pliku>
Na przykład:

Bash

make run INPUT=sinus
Czyszczenie
Aby wyczyścić miejsce robocze (usunąć skompilowane pliki), wykonaj:

Bash

make clean
Resetowanie
Aby przywrócić miejsce robocze do stanu początkowego (usunąć również wygenerowane dane), wykonaj:

Bash

make reset
