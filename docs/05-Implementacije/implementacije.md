# Linkovi za implementaciju

## GitHub repozitoriji (glavni)

Sljedeći repozitoriji sadrže postojeće implementacije simulacije nestišljivog fluida, uglavnom zasnovane na Stable Fluids algoritmu Jos Stama. Koristit ćemo ih kao reference za organizaciju simulacije, numeričke metode i optimizaciju performansi.

| # | Repozitorij | Jezik / tehnologija | Zašto je koristan |
|---|-------------|---------------------|-------------------|
| 1 | [torresds/fluids](https://github.com/torresds/fluids) | C++, OpenGL | 2D simulacija nestišljivog fluida zasnovana na Stable Fluids algoritmu. Korisna kao referenca za sekvencijalnu implementaciju i vizualizaciju kretanja fluida. |
| 2 | [woeishi/StableFluids](https://github.com/woeishi/StableFluids) | C++, SIMD (SSE/AVX) | Implementacija Stable Fluids algoritma s 2D i 3D verzijama i SIMD optimizacijama. Posebno korisna za analizu vektorizacije, organizacije podataka i korištenja Jacobijeve metode umjesto Gauss-Seidelove radi efikasnije SIMD obrade. |
| 3 | [taichi-dev/taichi – Stable Fluid](https://github.com/taichi-dev/taichi/blob/master/python/taichi/examples/graph/stable_fluid_graph.py) | Python, Taichi, Vulkan | 2D simulacija koja uključuje advekciju, računanje divergencije, Jacobijeve iteracije za pritisak i korekciju brzina. Korisna za razumijevanje koraka algoritma i vizualizaciju obojenog fluida. |
| 4 | [nornagon/stam-stable-fluids](https://github.com/nornagon/stam-stable-fluids) | JavaScript | Implementacija zasnovana na C kodu iz rada *Real-Time Fluid Dynamics for Games*. Korisna za razumijevanje jednostavnije strukture algoritma, dodavanja sila i izvršavanja vremenskih koraka simulacije. |

## Dodatne implementacije (vizualizacija)

- [PavelDoGreat/WebGL-Fluid-Simulation](https://github.com/PavelDoGreat/WebGL-Fluid-Simulation) – WebGL implementacija simulacije fluida u pregledniku. Korisna kao vizualna referenca za prikaz kretanja i miješanja obojene materije.

## Primjena u našem projektu

Planiramo implementirati 2D simulaciju nestišljivog fluida u C++-u, s automatskim pokretanjem fluida pomoću unaprijed definisanih sila.

Navedene implementacije poslužit će kao reference za razvoj sekvencijalnog algoritma, razumijevanje numeričkih postupaka i primjenu SIMD optimizacija. Nakon toga planiramo paralelizaciju odgovarajućih dijelova simulacije korištenjem OpenMP-a te poređenje performansi različitih verzija.