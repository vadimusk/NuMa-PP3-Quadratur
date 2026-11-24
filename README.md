# Programmierprojekt 3 – Numerische Integration (Kondensator-Energie)

Numerische Mathematik für Elektrotechnik, WS 2026/27. MATLAB (getestet mit R2026a, keine Toolboxen nötig).

## Ausführen

MATLAB im Projektordner starten und nacheinander ausführen:

```matlab
aufgabe1   % Aufgabe 1 + 2: Gauß und Newton-Cotes, Exaktheitsgrad, Kondition
aufgabe3   % Summierte Quadratur, Budget 48 Auswertungen
aufgabe4   % Kondensator: Quadratur- und Datenfehler, minimale Punktzahl
aufgabe5   % Experiment D: periodischer Integrand, Trapez gegen Gauß
```

Die Plots landen in `plots/`.

## Dateien

| Datei | Inhalt (Schnittstelle wie im Aufgabenblatt) |
|---|---|
| `gauss_formel.m` | `[x,w] = gauss_formel(N)` – Golub-Welsch über das Eigenwertproblem (Prop. 3.24) |
| `nc_formel.m` | `[x,w] = nc_formel(N)` – Newton-Cotes, Gewichte = ∫ℓᵢ, ℓᵢ baryzentrisch, Integral mit Gauß |
| `transformiere.m` | `[x,w] = transformiere(xh,wh,a,b)` – Lemma 3.3 |
| `summierte_formel.m` | `[x,w] = summierte_formel(xh,wh,a,b,M)` – M Teilintervalle, gemeinsame Randpunkte einmal |
| `quadratur.m` | `Q = quadratur(y,w)` |
| `exaktheitsgrad.m` | `m = exaktheitsgrad(xh,wh,tol)` – Monome xᵏ bis Abweichung > tol |
| `aufgabe1.m` … `aufgabe5.m` | Experimente und Plots |
| `speichere_plot.m` | Plot als PNG (immer helles Design) |

## Ergebnisse

- **Aufgabe 1/2:** Gauß: gemessener Grad = 2N−1, Σ|wᵢ| = 2 für alle N (alle wᵢ > 0). Newton-Cotes: Grad N (N ungerade) bzw. N−1 (N gerade), aber Σ|wᵢ| wächst exponentiell (N=20: 126, N=40: 1.6·10⁷), Vorzeichenwechsel in den Gewichten. Ab N≈19 (tol 1e-14) bzw. N≈34 (tol 1e-10) zeigen die Monome bei NC nur noch Rundungsverstärkung (ε·Σ|w|). Bei Gauß liegt der gemessene Grad für großes N *über* 2N−1: Der Fehler für x^(2N) ist ∝ 4⁻ᴺ und liegt unter der Rundungsgrenze, der Monomtest sieht ihn nicht.
- **Aufgabe 3:** f₁ (glatt): Gauß N≥4–5 erreicht Maschinengenauigkeit, NC braucht N≈8–10 und wird für N>16 durch Σ|w| schlechter. f₂ (schmaler Peak bei 0.3, Pole bei 0.3 ± 0.05i): Gauß mit mittlerem N (N=16, M=3: 9·10⁻⁸); sehr große N und sehr kleine N sind schlechter, NC wird bei großem N unbrauchbar. f₃ (Knick bei π/5): algebraische Konvergenz O(h²) unabhängig von N, Fehler ~10⁻⁴–10⁻⁵ bleibt; feinere Unterteilung (großes M) hilft mehr als hohe Ordnung.
- **Aufgabe 4:** W = 320·ln(1+144/64) = 377.1696 µJ. Mit 3 Gauß-Punkten: Quadraturfehler −0.48 µJ + Datenfehler 0.05 µJ, Worst-Case-Schranke für den Datenfehler 0.36 µJ → garantiert < 1 µJ (0.48 + 0.36 = 0.84). Newton-Cotes braucht 5, Simpson summiert 7 (garantiert), Trapez summiert 16 (garantiert, 12 real).
- **Aufgabe 5 (D):** Trapezregel mit N Punkten über eine Periode konvergiert geometrisch, Faktor ≈ e^(−1.62) ≈ 0.2 pro Punkt (aus den Polen von C(v(τ)), Streifenbreite s = 1.62; gemessen e^(−1.48 N) bei Phase τ₀ = 1). Fehler dauerhaft < 10⁻¹⁰ ab N = 21 (Trapez) gegen N = 32 (Gauß). Achtung: Startet die Abtastung bei τ = 0, sind wegen g(τ) = G(sin τ) die sin-Anteile der Fourier-Reihe unsichtbar, und für ungerades N fällt der Fehler wie e^(−2sN) (gemessen e^(−3.24 N), z. B. 3·10⁻¹¹ bei N = 11, aber 1.6·10⁻⁴ bei N = 12). Das ist ein Symmetrie-Artefakt, deshalb wird mit allgemeiner Phase bewertet.
