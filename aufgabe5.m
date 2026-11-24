%% Programmierprojekt 3 – Aufgabe 5: Experiment D (periodischer Integrand)
% Frage: Was passiert, wenn der Integrand periodisch ist? Trapezregel gegen Gauß.
%
% Beispiel aus der Aufgabe: Spannung v(t) = 6 V + 4 V sin(wt), Strom i = C(v) * v'.
% Mittelwert des Stromquadrats über eine Periode:
%     <i^2> = 1/T int_0^T i(t)^2 dt.
% Mit tau = w*t ist  i/w = C(v(tau)) * 4 cos(tau)  und  <i^2>/w^2 = 1/(2pi) int_0^{2pi} g(tau) dtau,
%     g(tau) = ( C(6 + 4 sin tau) * 4 cos tau )^2.
% Vermutung: Trapezregel (N Punkte, h = 2pi/N, Endpunkt nur einmal) konvergiert geometrisch,
% weil bei periodischen Funktionen die Randterme der Euler-Maclaurin-Formel wegfallen.

clear; close all; clc;

C = @(v) 10 ./ (1 + (v/8).^2);                       % uF
g = @(tau) (C(6 + 4*sin(tau)) .* 4 .* cos(tau)).^2;  % (uF*V)^2
T = 2*pi;

% Referenzwert: adaptive Integration, Kontrolle mit sehr feiner Trapezregel
Iref  = integral(g, 0, T, 'AbsTol', 1e-14, 'RelTol', 1e-13) / T;
Iref2 = mean(g((0:4095)' * T/4096));
fprintf('Referenz <i^2>/w^2 = %.12f  (Trapez N=4096: %.12f, Differenz %.1e)\n\n', Iref, Iref2, abs(Iref - Iref2));

% Achtung Symmetrie: g(tau) = G(sin tau), daher enthält die Fourier-Reihe nur cos(2m tau) und
% sin((2m+1) tau). Startet die Abtastung bei tau = 0, liegen alle Punkte in den Nullstellen der
% sin-Anteile -> diese sind unsichtbar, und für ungerades N ist die erste "aliasende" Frequenz 2N
% statt N (Fehler ~ exp(-2sN) statt exp(-sN)). Das ist ein Artefakt der Phase, kein Vorteil der
% Methode. Für die Bewertung wird deshalb eine allgemeine Phase tau0 verwendet.
tau0 = 1;

Ns = 1:60;
eT = nan(size(Ns));  eT0 = nan(size(Ns));  eG = nan(size(Ns));  eN = nan(size(Ns));
for N = Ns
    % Periodische Trapezregel: N Punkte tau_j = tau0 + j*T/N, alle Gewichte T/N
    tau = (0:N-1)' * T / N;
    eT(N)  = abs(mean(g(tau0 + tau)) - Iref);   % allgemeine Phase
    eT0(N) = abs(mean(g(tau)) - Iref);          % Start bei tau = 0 (symmetrische Abtastung)

    % Gauß-Legendre mit N Punkten auf [0, T]
    [xh, wh] = gauss_formel(N);
    [x, w]   = transformiere(xh, wh, 0, T);
    eG(N) = abs(quadratur(g(x), w) / T - Iref);

    % Newton-Cotes mit N Punkten (nur zum Vergleich, N >= 2, bis 40)
    if N >= 2 && N <= 40
        [xh, wh] = nc_formel(N);
        [x, w]   = transformiere(xh, wh, 0, T);
        eN(N) = abs(quadratur(g(x), w) / T - Iref);
    end
end

% Theorie: g ist analytisch im Streifen |Im tau| < s. Pole: 1 + (v/8)^2 = 0, also v = 6 + 4 sin(tau) = +-8i.
tau_pol = [asin((8i - 6)/4), asin((-8i - 6)/4)];
s = min(abs(imag(tau_pol)));
fprintf('Konvergenzstreifen: s = %.4f  ->  Trapezfehler ~ exp(-s*N): Faktor %.3f pro Punkt\n', s, exp(-s));

% Gemessene Rate (Fit für N >= 14, oberhalb des Rundungsbodens; davor noch nicht asymptotisch)
bereich = Ns(eT > 1e-12 & Ns >= 8);
p = polyfit(bereich, log(eT(bereich)), 1);
fprintf('gemessene Rate Trapez (tau0 = %g):         exp(%.3f*N)  (Faktor %.3f pro Punkt)\n', tau0, p(1), exp(p(1)));
bereich0 = Ns(eT0 > 1e-12 & Ns >= 5 & mod(Ns, 2) == 1);
p0 = polyfit(bereich0, log(eT0(bereich0)), 1);
fprintf('gemessene Rate Trapez (tau0 = 0, N ungerade): exp(%.3f*N)  (Theorie exp(-2sN) = exp(%.3f*N))\n', p0(1), -2*s);
fprintf('Trapez tau0 = 0, N gerade: Fehler(2N) = Fehler(N ungerade), z. B. e(10) = %.2e, e(5) = %.2e\n', eT0(10), eT0(5));
bereichG = Ns(eG > 1e-12 & Ns >= 14);
pG = polyfit(bereichG, log(eG(bereichG)), 1);
fprintf('gemessene Rate Gauß:   exp(%.3f*N)  (Faktor %.3f pro Punkt)\n\n', pG(1), exp(pG(1)));

fprintf('  N | Trapez     | Trapez t0=0 | Gauß       | Newton-Cotes\n');
for N = [4 5 8 9 11 12 16 20 24 28 32 40]
    fprintf('%3d | %.3e | %.3e   | %.3e | %.3e\n', N, eT(N), eT0(N), eG(N), eN(N));
end
% Punkte, ab denen der Fehler dauerhaft unter tol bleibt
ab = @(e, tolx) find(e >= tolx, 1, 'last') + 1;
for tolx = [1e-6 1e-10 1e-12]
    fprintf('Fehler dauerhaft < %.0e:  Trapez N = %d,  Gauß N = %d\n', tolx, ab(eT, tolx), ab(eG, tolx));
end

figure('Name', 'Aufgabe 5 (D)', 'Position', [100 100 1000 450]);
semilogy(Ns, max(eT, 1e-17), 'o-', 'LineWidth', 1.3); hold on;
semilogy(Ns, max(eT0, 1e-17), '.:', 'LineWidth', 1.0, 'MarkerSize', 10);
semilogy(Ns, max(eG, 1e-17), 's-', 'LineWidth', 1.3);
eNp = eN; eNp(~isnan(eN)) = max(eN(~isnan(eN)), 1e-17);
semilogy(Ns, eNp, '^-', 'LineWidth', 1.0);
semilogy(Ns, eT(20) * exp(-s*(Ns - 20)), 'k--');
grid on; ylim([1e-17 1e3]);
xlabel('N (Funktionsauswertungen)'); ylabel('Fehler von <i^2>/\omega^2  in (\muF V)^2');
legend(sprintf('Trapez, Phase \\tau_0 = %g', tau0), 'Trapez, Start bei \tau = 0', 'Gauß', 'Newton-Cotes', sprintf('exp(-%.2f N)', s), 'Location', 'southwest');
title('Periodischer Integrand g(\tau): Trapezregel gegen Gauß');
speichere_plot('a5_periodisch.png');
