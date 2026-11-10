# Linkovi za implementaciju

## GitHub repozitoriji

Sljedeći repozitoriji predstavljaju postojeće implementacije simulacije fluida zasnovane na PIC/FLIP metodama ili srodnim numeričkim pristupima. Odabrani su na osnovu relevantnosti za naš projekat, dostupnosti izvornog koda i zastupljenosti u open-source zajednici.

Koristit ćemo ih kao reference za organizaciju C++ koda, implementaciju numeričkih metoda, obradu čestica, simulacijsku mrežu i optimizaciju performansi.

| # | Repozitorij | Jezik / tehnologija | Zašto je koristan |
|---|---|---|---|
| 1 | [doyubkim/fluid-engine-dev](https://github.com/doyubkim/fluid-engine-dev) | C++, Python, OpenMP, TBB | Razvijen okvir za simulaciju fluida koji podržava 2D i 3D PIC, FLIP i APIC metode. Sadrži implementacije numeričkih solvera, prijenosa podataka između čestica i mreže te paralelne obrade. Posebno je koristan za našu C++ implementaciju i kasnije optimizacije. |
| 2 | [rlguy/Blender-FLIP-Fluids](https://github.com/rlguy/Blender-FLIP-Fluids) | C++, Python, Blender | Popularan projekat za simulaciju tečnosti zasnovan na FLIP metodi. Koristan je za proučavanje organizacije složenijeg simulacijskog sistema, rada s česticama i vizualizacije rezultata. Projekat je prvenstveno usmjeren na 3D simulaciju. |
| 3 | [thunil/mantaflow](https://github.com/thunil/mantaflow) | C++, Python | Okvir za istraživanje simulacije fluida s paralelizovanom C++ jezgrom. Koristan za razumijevanje organizacije numeričkih proračuna, upravljanja simulacijom i razdvajanja računarskog dijela od korisničkog interfejsa. |
| 4 | [davrempe/2d-fluid-sim](https://github.com/davrempe/2d-fluid-sim) | C++ | Implementacija dvodimenzionalne simulacije fluida zasnovane na PIC/FLIP pristupu. Korisna za proučavanje MAC mreže, prijenosa brzina između čestica i mreže, projekcije pritiska i organizacije sekvencijalne simulacije. |

## Dodatne implementacije i tutorijali

- [Matthias Müller – Ten Minute Physics, Tutorial 18](https://matthias-research.github.io/pages/tenMinutePhysics/index.html) – Praktični tutorijal za implementaciju FLIP simulacije vode. Koristan za razumijevanje osnovnih faza algoritma, uključujući prijenos brzina između čestica i mreže, proračun pritiska i ažuriranje položaja čestica.

- [boelukas/pbs-sandfluid](https://github.com/boelukas/pbs-sandfluid) – Studentski projekat koji implementira PIC i FLIP metode koristeći Taichi. Koristan za razumijevanje razlika između PIC i FLIP pristupa i poređenje rezultata simulacije.

## Primjena u našem projektu

U okviru projekta *Incompressible Fluid Simulation* implementiramo dvodimenzionalnu simulaciju nestišljivog fluida u C++-u korištenjem kombinovanog PIC/FLIP pristupa.

Simulacija se zasniva na česticama koje predstavljaju fluid i MAC mreži na kojoj se izvršavaju numerički proračuni.

Osnovni koraci simulacije uključuju:

1. Inicijalizaciju čestica i simulacijske mreže.
2. Prijenos brzina s čestica na mrežu (Particle-to-Grid).
3. Primjenu vanjskih sila, prvenstveno gravitacije.
4. Rješavanje jednačine pritiska i korekciju brzina.
5. Prijenos brzina s mreže nazad na čestice (Grid-to-Particle).
6. Ažuriranje položaja čestica i obradu sudara s granicama.

Navedene implementacije koristit ćemo kao reference za razumijevanje algoritma i provjeru vlastitog pristupa.

Prvo razvijamo sekvencijalnu C++ implementaciju, koja će predstavljati osnovu za poređenje performansi.

Nakon provjere numeričke ispravnosti planiramo implementirati SIMD vektorizaciju, OpenMP paralelizaciju i njihovu kombinaciju.

Cilj je analizirati vrijeme izvršavanja, ostvareno ubrzanje i efikasnost različitih implementacija pri jednakim početnim uslovima simulacije.