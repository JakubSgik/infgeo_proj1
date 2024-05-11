**README test** 

- udokumentowany, czyli opisany na githubie w pliku README.md, który ma zawierać następujące informacje:
    - do czego służy program i jaką funkcjonalność oferuje (transformacja XYZ -> BLH, transformacja BLH -> XYZ, jakie elipsidy są obsługiwane, ...)
    - jakie wymagania trzeba spełnić, by program działał na danym komputerze (np. trzeba mieć pythona w wesji takiej-a-takiej, 
      zainstalowaną bibliotekę taką-a-taką, ...)
    - dla jakiego systemu operacyjnego został napisany program 
    - jak go używać wraz z kilkoma przykładami wywołań obrazującymi jak z niego korzystać (w tym opis struktury danych wejściowych i wyjściowych) 
      oraz rezultatami tych wywołań (przykładowe wywołania powinny za input brać plik z przykładowymi danymi)
    - znane błędy i nietypowe zachowania programu, które nie zostały jeszcze naprawione
    
    
# p1_geo.py

p1_geo to nieduża biblioteka umożliwiająca podstawowe transformacje między układami współrzędnych używanymi w geodezji.

---

## Dostępne transformacje

W bibliotece zostały zaimplementowane jak dotąd metody pozwalające na następujące operacje:
- XYZ -> PLH
- PLH -> XYZ
- XYZ -> NEU
- PLH -> PL-1992
- PLH -> PL-2000

Ponadto, biblioteka zawiera wbudowane parametry wybranych elipsoid:
- WGS84     `wgs84`
- GRS'80    `grs80`
- Krasowski `krasowski`
- Mars      `mars`

```bash
pip install foobar
```

## Użycie


```bash
python p1_geo.py --xyz2plh wsp_inp.txt
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)