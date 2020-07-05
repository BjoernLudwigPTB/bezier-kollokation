import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.util.FastMath;

/**
 * Implementierung der in der Ausarbeitung beschriebenen Kollokationsmethode
 * unter Verwendung der Bernsteinbasis.
 */
public class BezierKollokation {
    /**
     * Die Genauigkeit, mit der gerechnet werden soll.
     */
    private static final int genauigkeit = 64;
    /**
     * Der verwendete Körper in dem gerechnet wird.
     */
    private final DfpField koerper = new DfpField(genauigkeit);
    /**
     * Die Anzahl der Kollokationspunkte je Gitterintervall.
     */
    private final int k;
    /**
     * Die Kollokationspunkte $\tau_1, ..., \tau_{l \cdot k}$.
     */
    private final Dfp[] tau;
    /**
     * Die Gitterknoten $s = \xi_0, \xi_1, ..., \xi_{l-1}, \xi_l = t$.
     */
    private final Dfp[] xi;
    /**
     * Die häufig auftretenden Faktoren $\mu_j^r := ((\tau_j-s)/(t - s))^r,$
     * $j = 1, ..., k, \quad r = 1, ..., k+1$. Es ist $mu[j-1][r-1]=\mu_j^r$.
     */
    private final Dfp[][] mu;
    /**
     * Die Näherungslösung $g$ von $-\varepsilon y'' - p y' + q y = f$.
     */
    private final BezierSplineFunction g;

    /**
     * Erzeugt eine Instanz des Näherungsverfahrens und berechnet die
     * Näherungslösung, wenn schon ein Gitter erzeugt ist.
     *
     * @param k       Anzahl der Kollokationspunkte je Gitterintervall.
     * @param xi      das Gitter auf dem die Näherungslösung berechnet werden soll.
     * @param epsilon singulärer Störungsparameter.
     * @param eta1    linker Randwert.
     * @param eta2    rechter Randwert.
     * @param p       Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param q       Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param f       in $-\varepsilon y'' - p y' + q y = f$.
     */
    public BezierKollokation(int k, Gitter<Dfp> xi, Dfp epsilon,
                             Dfp eta1, Dfp eta2, RealFieldPolynomialFunction p,
                             RealFieldUnivariateFunction<Dfp> q,
                             RealFieldUnivariateFunction<Dfp> f) {
        this.k = k;
        this.xi = xi.getXi();
        int l = this.xi.length - 1;
        tau = initialisiereTau(l);
        mu = initialisiereMus(l);
        Dfp[][] a = initialisiereA(l, epsilon, p, q);
        Dfp[] v = initialisiereV(l, eta1, eta2, f);
        Dfp[] loesung = new FieldBlockDecomposition(a, k, l)
                .solve(v);
        BezierFunction[] functions = new BezierFunction[l];
        for (int i = 0; i < l; i++) {
            functions[i] = new BezierFunction(koerper, ArrayUtils
                    .subarray(loesung, i * (k + 2), (i + 1) * (k + 2)),
                    getXi(i), getXi(i + 1));
        }
        g = new BezierSplineFunction(getXi(), functions);
    }

    /**
     * Erzeugt eine Instanz des Näherungsverfahrens und berechnet die
     * Näherungslösung auf einem uniformen Gitter.
     *
     * @param k    Anzahl der Kollokationspunkte je Gitterintervall.
     * @param l    die Anzahl der Gitterintervalle.
     * @param s    linkes Intervallende.
     * @param t    rechtes Intervallende.
     * @param eta1 linker Randwert.
     * @param eta2 rechter Randwert.
     * @param p    Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param q    Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param f    in $-\varepsilon y'' - p y' + q y = f$.
     */
    public BezierKollokation(int k, int l, Dfp s, Dfp t, Dfp eta1, Dfp eta2,
                             RealFieldUnivariateFunction<Dfp> p,
                             RealFieldUnivariateFunction<Dfp> q,
                             RealFieldUnivariateFunction<Dfp> f) {
        this.k = k;
        xi = initialisiereXi(l, s, t);
        tau = initialisiereTau(l);
        mu = initialisiereMus(l);
        Dfp[][] a = initialisiereA(l, koerper.getOne().negate(), p, q);
        Dfp[] v = initialisiereV(l, eta1, eta2, f);
        Dfp[] loesung = new FieldBlockDecomposition(a, k, l)
                .solve(v);
        BezierFunction[] functions = new BezierFunction[l];
        for (int i = 0; i < l; i++) {
            functions[i] = new BezierFunction(koerper, ArrayUtils
                    .subarray(loesung, i * (k + 2), (i + 1) * (k + 2)),
                    getXi(i), getXi(i + 1));
        }
        g = new BezierSplineFunction(getXi(), functions);
    }

    /**
     * Gibt die verwendete Genauigkeit zurück.
     *
     * @return die Anzahl der Nachkommastellen.
     */
    public static int getGenauigkeit() {
        return genauigkeit;
    }

    /**
     * Erzeugt für ein übergebenes $k$ und ein $l$ eine Instanz des
     * Kollokationsverfahrens für den klassischen Fall.
     *
     * @param k    die Anzahl der Kollokationspunkte je Gitterintervall.
     * @param l    Anzahl der Gitterintervalle.
     * @param mode Gibt an, welche Testdaten auf der Konsole ausgegeben werden
     *             sollen:
     *             $\verb!mode = 1!$: Gibt die Differenzen der Ableitungswerte an einigen
     *             Stellen des Intervalls $[0, 1]$ aus.
     *             $\verb!mode = 2!$: Gibt die Maplebefehle zum Plotten der Näherungslösung aus.
     *             $\verb!mode = 3!$: Gibt die Abweichungen von den Kollokations- und
     *             Randbedingungen aus.
     *             $\verb!mode = 4!$: Gibt die $\mu_j^i$ für alle zulässigen Parameterwerte aus.
     *             $\verb!mode = 5!$: Gibt alle $\tau_j$ aus.
     *             $\verb!mode = 6!$: Gibt alle Abweichungen der Ableitungen $i$-ter
     *             Ordnung, $i = 0, 1, 2$, der korrespondierenden exakten Lösung und der
     *             Näherungslösung aus.
     *             $\verb!mode = 7!$: Gibt alle $\xi_j$ aus.
     *             $\verb!mode = 8!$: Berechnet die Fehler und die experimentelle
     *             Konvergenzordnung für das klassische Beispiel mit $l = 2^7, \hdots, 2^{12}$.
     *             Ein Aufruf in einer Schleife ist nicht nötig.
     * @return Ausgabe gemäß $\verb!mode!$.
     */
    public static String berechneKlassisch(int k, int l, int mode) {
        StringBuilder ausgabe = new StringBuilder("k = " + k + ", l = " + l + "\n");
        /* Der Körper auf dem die Funktionen definiert sind. */
        final DfpField koerper = new DfpField(genauigkeit);
        /* Der linke Rand $s$ des Kollokationsintervalls $[s, t]$. */
        final Dfp s = koerper.getZero();
        /* Der rechte Rand $t$ des Kollokationsintervalls $[s, t]$. */
        final Dfp t = koerper.getOne();
        /* Der Randwert $y(s) = \eta_1$. */
        final Dfp eta1 = koerper.getZero();
        /* Der Randwert $y(t) = \eta_2$. */
        final Dfp eta2 = koerper.getZero();
        /* Die Koeffizientenfunktion $p$ in $-\varepsilon y'' - p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction(
                new Dfp[]{koerper.getZero()});
        /* Die Koeffizientenfunktion $q$ in $-\varepsilon y'' - p y' + q y = f$. */
        final RealFieldPolynomialFunction q = new RealFieldPolynomialFunction(
                new Dfp[]{koerper.getTwo().multiply(koerper.getTwo())
                        .negate()});
        /* Die Funktion $f$ in $-\varepsilon y'' - p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> f = new
                RealFieldUnivariateFunction<>() {

                    private final Dfp value = (koerper.getE().add(koerper.getE()
                            .reciprocal())).multiply(koerper.getTwo());

                    public Dfp value(Dfp x) {
                        return value;
                    }
                };
        BezierKollokation bkol = new BezierKollokation(k, l, s, t, eta1,
                eta2, p, q, f);
        BezierSplineFunction g = bkol.getG();
        final KlassischesBeispiel u;
        Dfp x, tMinusSDivN;
        int n;
        switch (mode) {
            case 1:
                n = 50;
                u = new KlassischesBeispiel(koerper);
                double xD;
                final double sD = s.toDouble();
                for (int i = 0; i <= n; i++) {
                    xD = sD + i * (t.toDouble() - sD) / n;
                    x = koerper.newDfp(xD);
                    for (int j : new int[]{0, 1, 2}) {
                        ausgabe.append("u^(").append(j).append(")(").append(x).append(") - g^(").append(j).append(")(").append(x).append(") = ").append(u.getAbleitung(xD, j) -
                                g.derivative(x, j).toDouble());
                        if (j < 2) ausgabe.append("\n");
                    }
                    if (i < n) ausgabe.append("\n");
                }
                break;
            case 2:
                n = 200;
                StringBuilder werte = new StringBuilder();
                StringBuilder stellen = new StringBuilder();
                tMinusSDivN = (t.subtract(s)).divide(n);
                for (int i = 0; i <= n; i++) {
                    x = s.add(tMinusSDivN.multiply(i));
                    werte.append(g.value(x).toDouble());
                    stellen.append(x.toDouble());
                    if (i != n) {
                        werte.append(",");
                        stellen.append(",");
                    }
                }
                ausgabe = new StringBuilder("u" + k + l + ":=pointplot([" + stellen + " ], [" + werte +
                        "], legend = \"u" + k + l + "\", color = \"" +
                        Linienfarbe.values()[(l + k - 1) % Linienfarbe.values().length]
                        + "\", connect = true):");
                break;
            case 3:
                ausgabe.append("g(").append(s).append(") = ").append(g.value(s)).append("\n");
                for (int j = 1; j <= l * k; j++) {
                    ausgabe.append("tau").append(j).append(" = ").append(bkol.getTau(j)).append(": g'' + p * g' + q * g - f = ").append(g.derivative(bkol.getTau(j), 2).subtract
                            (p.value(bkol.getTau(j)).multiply
                                    (g.derivative(bkol.getTau(j), 1))).add
                            (q.value(bkol.getTau(j)).multiply
                                    (g.derivative(bkol.getTau(j), 0))).subtract
                            (f.value(bkol.getTau(j)))).append("\n");
                }
                ausgabe.append("g(").append(t).append(")= ").append(g.value(t));
                break;
            case 4:
                for (int j = 1; j <= l * k; j++) {
                    for (int i = 0; i <= k + 1; i++) {
                        ausgabe.append("mu_").append(j).append("^").append(i).append(" = ").append(bkol.getMu(j / l, j % l, i, false)).append("\n");
                    }
                }
                break;
            case 5:
                for (int i = 0; i < l; i++) {
                    for (int j = 1; j <= k; j++) {
                        ausgabe.append("tau_").append(i * k + j).append(" = ").append(bkol.getTau(i,
                                j)).append("\n");
                    }
                }
                break;
            case 6:
                u = new KlassischesBeispiel(koerper);
                Dfp max = koerper.getZero();
                for (int j = 1; j <= l * k; j++) {
                    Dfp temp = g.value(bkol.getTau(j)).subtract(
                            u.value(bkol.getTau(j))).abs();
                    if (temp.greaterThan(max)) max = temp;
                }
                ausgabe.append("E_\\Delta^").append(k).append(Integer.toString(k).length() == 1 ?
                        " " : "").append(" = ").append(max).append("\n");
                max = koerper.getZero();
                for (int i = 0; i <= l; i++) {
                    Dfp temp = g.value(bkol.getXi(i)).subtract(
                            u.value(bkol.getXi(i))).abs();
                    if (temp.greaterThan(max)) max = temp;
                }
                ausgabe.append("E_\\xi^").append(l);
                ausgabe.append(" ".repeat(Math.max(0, Integer.toString(k).length() - Integer
                        .toString(l).length() - (Integer.toString(k)
                        .length() == 1 ? 0 : 1) + 4)));
                ausgabe.append(" = ").append(max).append("\n");
                break;
            case 7:
                Dfp[] tempXi = bkol.getXi();
                for (int j = 0; j <= l; j++) {
                    ausgabe.append("xi_").append(j + 1).append(" = ").append(tempXi[j]).append("\n");
                }
                break;
            case 8:
                ausgabe = new StringBuilder();
                int[] ls = new int[]{1, 2, 3, 4, 5, 6, 7, 8, 16, 32};
                /*
                 * max enthält die Fehler und experimentellen Konvergenzordnungen,
                 * wobei $\verb!max[0]! = E_\infty, \verb!max[1]! = \alpha_\infty,$ $\verb!max[2]! = E_\xi, \verb!max[3]! = \alpha_\xi,$
                 */
                Dfp[][] konvergenz = new Dfp[ls.length][4];
                for (int i = 0; i < ls.length; i++) {
                    konvergenz[i][0] = koerper.getZero();
                    konvergenz[i][2] = koerper.getZero();
                }
                /*
                 * Anzahl der pro Gitterintervall gleichverteilten Stellen, für
                 * welche die Näherung der Maximumnorm berechnet wird. Für
                 * Berechnungen die bis $k = 10, l = 32$ Änderungen der Norm
                 * kleiner $1^{-30}$ ergeben, sollte $n \geq 200$ sein.
                 */
                n = 50;
                u = new KlassischesBeispiel(koerper);
                for (int nu : new int[]{0, 1, 2}) {
                    System.out.println("k = " + k + ", nu = " + nu + ":\n");
                    for (int i = 0; i < ls.length; i++) {
                        int tempL = ls[i];
                        tMinusSDivN = (t.subtract(s)).divide(n * tempL);
                        bkol = new BezierKollokation(k, tempL, s,
                                t, eta1, eta2, p, q, f);
                        g = bkol.getG();
                        for (int j = 0; j <= n * tempL; j++) {
                            x = s.add(tMinusSDivN.multiply(j));
                            Dfp temp = g.derivative(x, nu).subtract(
                                    u.getAbleitung(x, nu)).abs();
                            if (temp.greaterThan(konvergenz[i][0])) {
                                konvergenz[i][0] = temp;
                            }
                        }
                        System.out.print("l = ");
                        int i1 = 2 - Integer.toString(ls[i]).length();
                        for (int m = 0; m < i1;
                             m++) {
                            System.out.print(" ");
                        }
                        System.out.println(ls[i] + ": E_\\infty = " +
                                konvergenz[i][0]);
                        if (i > 0) {
                            konvergenz[i][1] = konvergenz[i][0].divide(
                                    konvergenz[i - 1][0]).log().divide(
                                    FastMath.log((double) ls[i - 1] / ls[i]));
                            for (int m = 0; m < i1; m++) {
                                System.out.print(" ");
                            }
                            System.out.println("       \\alpha_" + ls[i] + " = " +
                                    konvergenz[i][1]);
                        }
                        for (int j = 0; j <= ls[i]; j++) {
                            Dfp temp = g.derivative(bkol.getXi(j), nu).subtract(
                                    u.getAbleitung(bkol.getXi(j), nu)).abs();
                            if (temp.greaterThan(konvergenz[i][2])) {
                                konvergenz[i][2] = temp;
                            }
                        }
                        System.out.println("           E_\\xi = " +
                                konvergenz[i][2]);
                        if (i > 0) {
                            konvergenz[i][3] = konvergenz[i][2].divide(
                                    konvergenz[i - 1][2]).log().divide(FastMath
                                    .log((double) ls[i - 1] / ls[i]));
                            for (int m = 0; m < i1;
                                 m++) {
                                System.out.print(" ");
                            }
                            System.out.println("        \\beta_" + ls[i] + " = " +
                                    konvergenz[i][3]);
                        }
                        System.out.println();
                    }
                }
                break;
            default:
                throw new IllegalStateException("Unexpected value: " + mode);
        }
        return ausgabe.toString();
    }

    /**
     * Berechnen verschiedener Daten in Abhängigkeit des angegebenen Modus
     * und Ausgabe auf der Konsole.
     */
    public static void main(String[] args) {
        int[] n = new int[]{1, 2, 3, 4, 5, 6, 7, 8, 16, 32};
        for (int k : new int[]{1, 2, 4}) {
            for (int l : n) {
                System.out.println(berechneKlassisch(k, l, 6));
            }
        }
    }

    /**
     * Berechnet alle $\mu_j^r := \mu(\tau_j)^r$ und $(1 - \mu_j)^r := (1 - \mu(\tau_j))^r$ und speichert diese.
     * Vor dem Aufruf dieser Prozedur muss das Feld $\verb!Dfp[] tau!$ mit den
     * Kollokationspunkten initialisiert sein. Die Berechnung der
     * Funktionswerte ist dabei auf ein Minimum beschränkt unter Ausnutzung
     * von $\mu_j = 1 - \mu_{k-j-1}$ und äquivalent $1 - \mu_j = \mu_{k-j-1}$.
     *
     * @param l die Anzahl der Gitterintervalle.
     * @return ein Feld $\verb!Dfp[] mu!$ mit $\verb!mu[j][r]! = \mu_j^r$.
     */
    private Dfp[][] initialisiereMus(int l) {
        Dfp[][] tempMu = new Dfp[l * k][k + 1];
        for (int i = 0; i < l; i++) {
            BezierFunction tempG = new BezierFunction(koerper,
                    new Dfp[]{koerper.getZero(), koerper.getZero()},
                    getXi(i), getXi(i + 1));
            int index1 = i * k;
            for (int j = 0; j < k; j++) {
                int index2 = index1 + j;
                if (j > k - j)
                    tempMu[index2][0] = koerper.getOne()
                            .subtract(tempMu[index1 + k - j - 1][0]);
                else
                    tempMu[index2][0] = tempG.getMu(getTau(index2 + 1));
                for (int r = 1; r <= k; r++) {
                    tempMu[index2][r] = tempMu[index2][r - 1]
                            .multiply(tempMu[index2][0]);
                }
            }
        }
        return tempMu;
    }

    /**
     * Berechnet die streng monoton steigende Folge $\tau_j, j = 1, ..., l \cdot k$.
     */
    private Dfp[] initialisiereTau(int l) {
        Dfp[] temptau = new Dfp[l * k];
        GaussLegendrePunkte rhos = new GaussLegendrePunkte(k, koerper);
        for (int i = 0; i < l; i++) {
            Dfp tempPlus = getXi(i).add(getXi(i + 1)), tempMinus = getXi(i + 1)
                    .subtract(getXi(i));
            for (int j = 0; j < k; j++) {
                temptau[(i * k) + j] = tempPlus.add(tempMinus.multiply(rhos
                        .getRho(j))).divide(koerper.getTwo());
            }
        }
        return temptau;
    }

    /**
     * Berechnet eine im Intervall $[s, t]$ uniforme Gitterknotenfolge
     * $\xi_i, i = 0, ..., l$ mit $\xi_0 = s,\xi_l = t$.
     *
     * @param s das linke Intervallende von $[s, t]$.
     * @param t das rechte Intervallende von $[s, t]$.
     * @param l die Anzahl der Gitterintervalle.
     */
    private Dfp[] initialisiereXi(int l, Dfp s, Dfp t) {
        return new UniformesGitter(l, s, t).getXi();
    }

    /**
     * Berechnet die Einträge der Koeffizientenmatrix $A$ des zu lösenden
     * Gleichungssystems.
     *
     * @param l die Anzahl der Gitterintervalle.
     * @param p in $-\varepsilon y'' - p y' + q y = f$.
     * @param q in $-\varepsilon y'' - p y' + q y = f$.
     * @return $A$.
     */
    private Dfp[][] initialisiereA(int l, Dfp epsilon,
                                   RealFieldUnivariateFunction<Dfp> p,
                                   RealFieldUnivariateFunction<Dfp> q) {
        /*
         * Für $k = 1$ entsteht bei der blockweisen Erstellung die
         * Besonderheit, dass die Blöcke der Länge $k + 2 = 3$ nicht mehr
         * die Stetigkeitsbedingungen aufnehmen können, da sich diese
         * jeweils auf vier Bézierpunkte beziehen. Deshalb wird in diesem
         * Fall die Matrix um eine Spalte erweitert. Eine detaillierte
         * Beschreibung der Blockstruktur von $A$ ist in
         * $\verb!FieldBlockDecomposition!$ zu finden.
         */
        Dfp[][] tempA = new Dfp[l * (k + 2)][(k == 1 && l > 1) ? k + 3 : k + 2];
        /* Variablen zum Zwischenspeichern wiederholt benötigter
         * Zwischenergebnisse und der auftretenden Binomialkoeffizienten.
         */
        Dfp deltaXi = koerper.getOne();
        Binomialkoeffizient binomM = new Binomialkoeffizient(k - 1),
                binomK = new Binomialkoeffizient(k),
                binomP = new Binomialkoeffizient(k + 1);
        /*
         * Befüllt die erste und letzte Zeile der Matrix $A$ mit den
         * Randbedingungen.
         */
        tempA[0][0] = koerper.getOne();
        for (int i = 1; i < tempA[0].length; i++) {
            tempA[0][i] = koerper.getZero();
        }
        for (int i = 0; i < tempA[0].length - 1; i++) {
            tempA[tempA.length - 1][i] = koerper.getZero();
        }
        tempA[tempA.length - 1][tempA[0].length - 1] = koerper.getOne();
        /* Erzeugen der Blockstruktur durch verschachtelte Schleifen. */
        for (int i = 0; i < l; i++) {
            {
                /* Wiederholt benötigte Ausdrücke. */
                Dfp deltaXiMinus = deltaXi;
                deltaXi = getXi(i + 1).subtract(getXi(i));
                /*
                 * Befüllt die $2$ Zeilen der Matrix $A$, welche aus den
                 * Stetigkeitsbedingungen für den Übergang vom $i-1$-ten zum
                 * $i$-ten Gitterintervall hervorgehen, $i = 1, \hdots, l-1$.
                 */
                if (i > 0) {
                    /* Die Stetigkeitsbedingung bezüglich der $C^1$-Stetigkeit. */
                    tempA[i * (k + 2) - 1][0] = deltaXi;
                    tempA[i * (k + 2) - 1][1] = deltaXiMinus.add(deltaXi).negate();
                    tempA[i * (k + 2) - 1][2] = koerper.getZero();
                    tempA[i * (k + 2) - 1][3] = deltaXiMinus;
                    /* Die Stetigkeitsbedingung bezüglich der $C^0$-Stetigkeit. */
                    tempA[i * (k + 2)][0] = koerper.getZero();
                    tempA[i * (k + 2)][1] = koerper.getOne();
                    tempA[i * (k + 2)][2] = koerper.getOne().negate();
                    tempA[i * (k + 2)][3] = koerper.getZero();
                    for (int j = 4; j < tempA[0].length; j++) {
                        tempA[i * (k + 2) - 1][j] = koerper.getZero();
                        tempA[i * (k + 2)][j] = koerper.getZero();
                    }
                }
            }
            /* Wiederholt benötigte Ausdrücke. */
            Dfp kPlusDivDeltaXi = deltaXi.reciprocal().multiply(k + 1),
                    epsilonKPlusKDivDeltaXiSqr = kPlusDivDeltaXi
                            .multiply(k).multiply(epsilon).divide(deltaXi),
                    deltaXiSqr = deltaXi.pow(2);
            /*
             * Befüllt die $k$ Zeilen der Matrix $A$, welche aus den
             * Kollokationsbedingungen des $i$-ten Gitterintervalls hervorgehen.
             */
            for (int j = 1; j <= k; j++) {
                /*
                 * Berechnung der Funktionswerte $p(\tau_{ik+j})$ und $q(\tau_{ik+j})$ für den
                 * $j$-ten Kollokationspunkt des $i$-ten Gitterintervalls.
                 */
                Dfp pJ = p.value(getTau(i, j)), qJ = q.value(getTau(i, j));
                /*
                 * Berechnung des Summanden bezüglich $b_{i0}$ der $j$-ten
                 * Kollokationsbedingung.
                 */
                Dfp[] dfps = tempA[i * (k + 2) + j];
                dfps[0] = getMu(i, j, k - 1, true)
                        .multiply(
                                pJ.multiply(kPlusDivDeltaXi)
                                        .multiply(
                                                getMu(i, j, 1, true)
                                        )
                                        .add(qJ.multiply(
                                                getMu(i, j, 2, true)
                                        ))
                                        .subtract(epsilonKPlusKDivDeltaXiSqr)
                        );
                /*
                 * Berechnung des Summanden bezüglich $b_{i1}$ der $j$-ten
                 * Kollokationsbedingung.
                 */
                dfps[1] = getMu(i, j, k - 2, true)
                        .multiply(k + 1).multiply(
                                epsilon.multiply(k)
                                        .divide(deltaXiSqr).multiply(
                                        koerper.getTwo()
                                                .subtract(getMu(i, j, 1, false)
                                                        .multiply(k + 1))
                                )
                                        .subtract(
                                                pJ.divide(deltaXi)
                                                        .multiply(
                                                                koerper.getOne()
                                                                        .subtract(
                                                                                getMu(i, j, 1, false)
                                                                                        .multiply(k + 1)
                                                                        )
                                                        )
                                                        .multiply(getMu(i, j, 1, true))
                                        )
                                        .add(
                                                qJ.multiply(getMu(i, j, 2, true))
                                                        .multiply(getMu(i, j, 1, false))
                                        )
                        );
                /*
                 * Berechnung des Summanden bezüglich $b_{i\kappa}, \kappa = 2, ..., k-1$ der $j$-ten
                 * Kollokationsbedingung.
                 */
                for (int kappa = 2; kappa < k; kappa++) {
                    dfps[kappa] =
                            epsilonKPlusKDivDeltaXiSqr
                                    .multiply(
                                            koerper.getTwo().multiply(
                                                    binomM.getUeber(kappa - 1))
                                                    .multiply(getMu(i, j, 1, true)
                                                            .multiply(getMu(i, j, 1, false)))
                                                    .subtract(getMu(i, j, 2, true)
                                                            .multiply(binomM
                                                                    .getUeber(kappa - 2)))
                                                    .subtract(getMu(i, j, 2, false).multiply(
                                                            binomM.getUeber(kappa)))
                                    )
                                    .multiply(getMu(i, j, k - 1 - kappa, true))
                                    .multiply(getMu(i, j, kappa - 2, false))

                                    .subtract(pJ.multiply(kPlusDivDeltaXi)
                                            .multiply(
                                                    getMu(i, j, 1, false)
                                                            .multiply(binomP.getUeber(kappa))
                                                            .negate()
                                                            .add(binomK.getUeber(kappa - 1))
                                            )
                                            .multiply(getMu(i, j, k - kappa, true))
                                            .multiply(getMu(i, j, kappa - 1, false)))

                                    .add(qJ.multiply(binomP.getUeber(kappa))
                                            .multiply(getMu(i, j, k + 1 - kappa, true))
                                            .multiply(getMu(i, j, kappa, false)));
                }
                /*
                 * Berechnung des Summanden bezüglich $b_{ik}$ der $j$-ten
                 * Kollokationsbedingung.
                 */
                dfps[k] = getMu(i, j, k - 2, false)
                        .multiply(k + 1).multiply(
                                epsilon.multiply(k)
                                        .divide(deltaXiSqr)
                                        .multiply(
                                                getMu(i, j, 1, false).multiply(k + 1)
                                                        .subtract(k).add(koerper.getOne())
                                        )
                                        .subtract(
                                                pJ.divide(deltaXi)
                                                        .multiply(
                                                                getMu(i, j, 1, false)
                                                                        .multiply(k + 1).negate()
                                                                        .add(k)
                                                        )
                                                        .multiply(getMu(i, j, 1, false))
                                        )
                                        .add(
                                                qJ.multiply(getMu(i, j, 2, false))
                                                        .multiply(getMu(i, j, 1, true))
                                        )
                        );
                /*
                 * Berechnung des Summanden bezüglich $b_{i,k+1}$ der $j$-ten
                 * Kollokationsbedingung.
                 */
                dfps[k + 1] = getMu(i, j, k - 1, false)
                        .multiply(
                                getMu(i, j, 1, false).multiply(
                                        qJ
                                                .multiply(
                                                        getMu(i, j, 1, false)
                                                )
                                                .subtract(
                                                        pJ.multiply(kPlusDivDeltaXi)
                                                )
                                )
                                        .subtract(epsilonKPlusKDivDeltaXiSqr)
                        );
                if (k == 1 && l > 1) dfps[k + 2] = koerper.getZero();
            }
        }
        return tempA;
    }

    /**
     * Berechnet die Einträge der rechten Seite $v$ des zu lösenden
     * Gleichungssystems.
     *
     * @param l    die Anzahl der Gitterintervalle.
     * @param eta1 Randwert $y(s) = \eta_1$.
     * @param eta2 Randwert $y(t) = \eta_2$.
     * @param f    in $-\varepsilon y'' + p y' + q y = f$.
     * @return $v$.
     */
    private Dfp[] initialisiereV(int l, Dfp eta1,
                                 Dfp eta2, RealFieldUnivariateFunction<Dfp> f) {
        Dfp[] tempV = new Dfp[l * (k + 2)];
        tempV[0] = eta1;
        for (int i = 0; i < l; i++) {
            int index = i * (k + 2);
            if (i > 0) {
                tempV[index - 1] = koerper.getZero();
                tempV[index] = koerper.getZero();
            }
            for (int j = 1; j <= k; j++) {
                tempV[index + j] = f.value(getTau(i, j));
            }
        }
        tempV[l * (k + 2) - 1] = eta2;
        return tempV;
    }

    /**
     * Gibt $\tau_j, j = 1, ..., l \cdot k$ zurück.
     *
     * @param j für das $\tau_j, j = 1, ..., l \cdot k$ zurückgegeben werden soll.
     * @return $\tau_j$
     */
    public Dfp getTau(int j) {
        return tau[j - 1];
    }

    /**
     * Gibt $\tau_j, j = 1, ..., k$ aus dem $i$-ten Gitterintervall zurück.
     *
     * @param i für das $i$-te Gitterintervall.
     * @param j für das $\tau_j, j = 1, ..., k$ zurückgegeben werden soll.
     * @return $\tau_{ik+j}$
     */
    public Dfp getTau(int i, int j) {
        return getTau(i * k + j);
    }

    /**
     * Gibt $\mu_j^{\operatorname{exponent}}$ oder $(1 - \mu_j)^{\operatorname{exponent}}$ zurück.
     *
     * @param i        das $i$-te Gitterintervall $i = 0, \hdots, l - 1$.
     * @param j        aus $1, ..., k$ zur Auswahl von $\mu_j$ oder $(1 - \mu_j)$.
     * @param exponent ${\operatorname{exponent}} \in \{0, ..., k+1\}$, für den $\mu_j^{\operatorname{exponent}}$ oder $(1 - \mu_j)^{\operatorname{exponent}}$.
     * @param invers   $\verb!true!$, falls $(1 - \mu_j)^{\operatorname{exponent}}$ und $\verb!false!$, falls $\mu_j^{\operatorname{exponent}}$
     *                 zurückgegeben werden soll.
     * @return $\mu_j^{\operatorname{exponent}}$ oder $(1 - \mu_j)^{\operatorname{exponent}}$
     * @throws ArrayIndexOutOfBoundsException falls ${\operatorname{exponent}} < 0$ oder ${\operatorname{exponent}} > k+1$.
     */
    private Dfp getMu(int i, int j, int exponent, boolean invers)
            throws ArrayIndexOutOfBoundsException {
        if (exponent < -1 || exponent > k + 1 || i < 0
                || i > xi.length + 1)
            throw new ArrayIndexOutOfBoundsException();
        if (exponent == -1)
            if (invers)
                return mu[(i + 1) * k - j][0].reciprocal();
            else
                return mu[i * k + j - 1][0].reciprocal();
        else if (exponent == 0)
            return koerper.getOne();
        else if (invers) {
            return mu[(i + 1) * k - j][exponent - 1];
        } else {
            return mu[i * k + j - 1][exponent - 1];
        }
    }

    /**
     * Gibt die bestimmte Näherungslösung zurück.
     *
     * @return g
     */
    public BezierSplineFunction getG() {
        return g;
    }

    /**
     * Gibt eine Kopie des Feldes $\verb!Dfp[] xi!$ der Gitterknoten
     * im Intervall $[s, t]$ zurück.
     *
     * @return eine Kopie von $\verb!Dfp[] xi!$
     */
    public Dfp[] getXi() {
        Dfp[] out = new Dfp[xi.length];
        System.arraycopy(xi, 0, out, 0, xi.length);
        return out;
    }

    /**
     * Gibt $\xi_i, i = 0, ..., l$ zurück.
     *
     * @param i für das $\xi_i, i = 0, ..., l$ zurückgegeben werden soll.
     * @return $\xi_i$
     */
    public Dfp getXi(int i) {
        return xi[i];
    }

    /**
     * Auszug aus den möglichen Farbnamen in Maple.
     */
    enum Linienfarbe {
        Brown, Crimson, Chocolate, Orange, SkyBlue, Magenta, Gold
    }
}
