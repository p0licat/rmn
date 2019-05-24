# rmn
Proiect in C++ pentru lucrul cu numere mari si de precizie arbitrara.

## Configurare

Proiectul necesita librariile GMP ( Gnu Multiple Precision Arithmetic ) si MPFR ( Multiple Precision Arithmetic with Correct Rounding ). Aceste librarii sunt Open-Source sub licenta GNU.

Pot fi gasite aici: [GMP](https://gmplib.org/), [MPFR](http://www.mpfr.org/).

Ele necesita compilate local si instructiunile pentru fiecare le gasiti pe site-urile respective. Trebuie instalate in ordinea:
```
1. GMPlib
2. MPFR
```

Proiectul de asemenea foloseste un wrapper pentru clasa MPFR scris de [Pavel Holoborodko](http://www.holoborodko.com/pavel/mpfr/), Open-Source sub licenta GNU, care este inclus si nu necesita instalare.

## Utilizare

Momentan programul are implementata o clasa functionala GmpMatrix, care aloca dinamic in memorie o matrice de valori "mpreal". Fiecare element al acestei matrici are o precizie stabilita de utilizator, si limitarile sunt doar memoria disponibila a sistemului pe care ruleaza.

Compilarea si executia proiectului in aceasta versiune lanseaza testele.









