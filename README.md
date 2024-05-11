# p1_geo.py

p1_geo.py to nieduża biblioteka napisana w języku Python, umożliwiająca podstawowe transformacje między układami współrzędnych używanymi w geodezji.



## Dostępne transformacje

W bibliotece zostały zaimplementowane jak dotąd metody pozwalające na następujące operacje:
- XYZ -> PLH `--xyz2plh`
- PLH -> XYZ `--plh2xyz`
- XYZ -> NEU `--xyz2neu`
- PLH -> PL-1992 `--plh1992`
- PLH -> PL-2000 `--plh2000`
- XYZ -> PL-1992 `--xyz1992`
- XYZ -> PL-2000 `--xyz2000`
---
Ponadto, biblioteka zawiera wbudowane parametry wybranych elipsoid:
- WGS84     `wgs84`
- GRS'80    `grs80`
- Krasowskiego `krasowski`
- Mars      `mars`

Domyślną elipsoidą jest **WGS84**. Aby zastosować model wybranej innej elipsoidy w obliczeniach, należy podać przy wywoływaniu programu flagę `--model [nazwa_modelu]` z nazwą wybranej elipsoidy, np.
```
--model krasowski
```
---
Program obsługuje również flagę `--naglowek [dlugosc_naglowka]`, która zapewnia poprawne działanie poprzez pominięcie podanej liczby początkowych wierszy przy odczytywaniu pliku z danymi wsadowymi. Domyślnie pomijany jest **1** wiersz, co również można skonfigurować przy wywołaniu:
```
--naglowek 5
```

## Przykłady użycia
Program uruchamia się w wierszu poleceń, podając każdorazowo nazwę biblioteki `p1_geo.py`, typ transformacji do wykonania oraz nazwę pliku z danymi wejściowymi (wskazanie pliku musi odbyć się zawsze na samym końcu polecenia).
```bash
python p1_geo.py [transformacja] [sciezka_do_pliku_z_danymi]
```
```bash
python p1_geo.py --xyz2plh wsp_inp.txt
```
Podanie długości nagłówka lub modelu elipsoidy jest opcjonalne. Nie ma znaczenia kolejność podawania flag - program automatycznie rozpozna odpowiednią flagę.

### Przykładowe polecenia
```bash
python p1_geo.py --xyz2plh --naglowek 4 --model grs80 wsp_inp.txt
```

```bash
python p1_geo.py --naglowek 2 --xyz1992 --model krasowski wsp_inp.txt
```
**Uwaga** - w przypadku transformacji XYZ -> NEU, współrzędne punktu środkowego powinny zostać podane w kolejności `X0 Y0 Z0` bezpośrednio przed podaniem pliku z danymi.
```bash
python p1_geo.py --model wgs84 --xyz2neu 3664940.500 1409153.590 5009571.170 wsp_inp.txt
```
### Struktura danych
W pliku z danymi wejściowymi może znajdować się nagłówek, którego liczba wierszy objętości powinna zostać wskazana przy wywołaniu dla prawidłowego działania programu.

Plik powinien być zapisany w formacie `.txt`.

Wewnątrz powinny znaleźć się dane liczbowe w kolejności wskazanej przez nazwę wybranej transformacji i rozdzielone przecinkami. Czyli np. dla `--xyz2plh` dane powinny być w zapisane w kolejności `X, Y, Z`. Analogicznie, dla `--plh2xyz` będzie to `phi, lam, h`.

Ta sama zasada znajduje zastosowanie przy generowaniu pliku z wynikami transformacji.


## Wymagania systemowe
Biblioteka została napisana dla systemu operacyjnego **Windows 11**.

Wersja Pythona: **3.11.7**

Niezbędne biblioteki:
- NumPy (wersja 1.26.4)

## Znane błędy
- brak możliwości podania współrzędnych `X0 Y0 Z0` dla transformacji do NEU w dowolnym miejscu 