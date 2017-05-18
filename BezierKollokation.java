import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.*;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.util.FastMath;

/**
 * Implementierung der in der Ausarbeitung beschriebenen Kollokationsmethode
 * unter Verwendung der Bernsteinbasis.
 */
public class BezierKollokation {
    /**
     * Die Genauigkeit, mit der gerechnet werden soll. Ab 18
     * Kollokationspunkten bringt eine Genauigkeit bis 16 Stellen immer das
     * gleiche Ergebnis.
     */
    private static final int genauigkeit = 17; 
    /**
     * Der verwendete Körper mit der Anzahl der Stellen, mit denen gerechnet
     * wird.
     */
    private DfpField koerper = new DfpField(genauigkeit);
    /** Die Randwerte $y(s) = eta1$ und $y(t) = eta2$. */
    private final Dfp eta1, eta2;
    /** Die Anzahl $k$ von Kollokationspunkten. */
    private final int k;
    /** Die Anzahl $l$ von Gitterknoten. */
    private final int l;
    /** Die Kollokationspunkte $\tau_1, ..., \tau_k$. */
    private final Dfp[] tau;
    /** Die Gitterknoten $s = \xi_1, \xi2, ..., \xi_l, \xi_{l+1} = t$. */
    private final Dfp[] xi;
    /**
     * Die häufig auftretenden Faktoren $\mu_j^i := ((\tau_j-s)/(t - s))^i,$
     * $j = 1, ..., k, \quad i = 1, ..., k+1$. Es ist $mu[j-1][i-1]=\mu_j^i$.
     */
    private final Dfp[][] mu;
    /** Koeffizientenmatrix $A$ des linearen Gleichungssystems. */
    private final SparseFieldMatrix<Dfp> A;
    /** rechte Seite $b$ des linearen Gleichungssystems. */
    private final SparseFieldVector<Dfp> v;
    /**
     * Die an der Differentialgleichung $y'' + p y' + q y = f$
     * beteiligten Funktionen.
     */
    private RealFieldUnivariateFunction<Dfp> p, q, f;
    /** Die Näherungslösung $g$ von $y'' + p y' + q y = f$. */
    private BezierDarstellung g;    
    
    /**
     * Erzeugt eine Instanz des Näherungsverfahrens und stößt die Berechnung
     * der Näherungslösung an.
     * @param k Anzahl der zu verwendenden Kollokationspunkte.
     * @param l Anzahl der zu verwendenden Gitterknoten $- 1$.
     * @param s linkes Intervallende.
     * @param t rechtes Intervallende.
     * @param eta1 linker Randwert.
     * @param eta2 rechter Randwert.
     * @param p Koeffizientenfunktion in $y'' + p y' + q y = f$.
     * @param q Koeffizientenfunktion in $y'' + p y' + q y = f$.
     * @param f rechte Seite der Differentialgleichung
     * $y'' + p y' + q y = f$.
     */
    public BezierKollokation(int k, int l, Dfp s, Dfp t, Dfp eta1, Dfp eta2,
            RealFieldUnivariateFunction<Dfp> p,
            RealFieldUnivariateFunction<Dfp> q,
            RealFieldUnivariateFunction<Dfp> f) {
        this.k = k;
        this.l = l;
        this.eta1 = eta1;
        this.eta2 = eta2;
        this.p = p;
        this.q = q;
        this.f = f;
        xi = initialisiereXi(s, t);
        tau = initialisiereTau();
        mu = initialisiereMus();
        A = initialisiereA();
        v = initialisiereV();
        FieldVector<Dfp> loesung = new FieldLUDecomposition<Dfp>(A)
                .getSolver().solve(v);
        g = new BezierDarstellung(koerper, loesung.toArray(), getS(), getT());
    }
    
    /**
     * Berechnet alle $\mu_j^i := \mu(\tau_j)^i$ und
     * $(1 - \mu_j)^i := (1 - \mu(\tau_j))^i$ und speichert diese. Vor
     * dem Aufruf dieser Prozedur muss das Feld {@code Dfp[] tau} mit den
     * Kollokationspunkten initialisiert sein. Die Berechnung der
     * Funktionswerte ist dabei auf ein Minimum beschränkt unter Ausnutzung
     * von <p>
     * $\mu_j = 1 - \mu_(k-j-1)$,
     * </p>
     * und äquivalent <p>
     * $1 - \mu_j = \mu_(k-j-1)$.</p>
     */
    private Dfp[][] initialisiereMus () {
        Dfp[][] tempmu = new Dfp[k][k+1];
        for (int j = 0; j < k; j++) {
            if (j > k - j)
                tempmu[j][0] = koerper.getOne()
                .subtract(tempmu[k - j - 1][0]);
            else
                tempmu[j][0] = BezierDarstellung
                .mu(getTau(j+1), getS(), getT());
            for (int i = 1; i <= k; i++) {
                tempmu[j][i] = tempmu[j][i-1].multiply(tempmu[j][0]);
            }
        }
        return tempmu;
    }
    
    /**
     * Berechnet die streng monoton steigende Folge
     * $\tau_j, j = 1, ..., k$ im Intervall $[s, t]$.
     */
    private Dfp[] initialisiereTau () {
        Dfp[] temptau = new Dfp[k];
        GaussLegendrePunkte rhos = new GaussLegendrePunkte(k, koerper);
        for (int j = 0; j < k; j++) {
            temptau[j] = (getS().add(getT()).add((getT().subtract(getS()))
                    .multiply(rhos.getRho(j)))).divide(koerper.getTwo());
        }
        return temptau;
    }
    
    /**
     * Berechnet die streng monoton steigende Folge
     * $\xi_j, j = 1, ..., l+1$ im Intervall $[s, t]$.
     */
    private Dfp[] initialisiereXi (Dfp s, Dfp t) {
        final Dfp c = t.subtract(s).divide(l);
        Dfp[] tempXi = new Dfp[l+1];
        Dfp temp;
        for (int j = 0; j <= FastMath.floor(l/2); j++) {
            temp = c.multiply(j);
            tempXi[j] = s.add(temp);
            tempXi[l - j] = t.subtract(temp);
        }
        return tempXi;
    }
    
    /**
     * Berechnet die Einträge der Koeffizientenmatrix $A$ des zu lösenden
     * Gleichungssystems.
     * @return $A$.
     */
    private SparseFieldMatrix<Dfp> initialisiereA () {
        /** Die Binomialkoeffizienten, deren Werte häufig gebraucht werden. */
        Binomialkoeffizient KMinusEins = new Binomialkoeffizient(k - 1);
        Binomialkoeffizient K = new Binomialkoeffizient(k);
        Binomialkoeffizient KPlusEins = new Binomialkoeffizient(k + 1);
        /** Häufig verwendete Berechnungsschritte und Funktionswerte. */
        Dfp TMinusS = getT().subtract(getS()),
                kPlusDivTMinusS = TMinusS.reciprocal().multiply(k + 1),
                kPlusKDivTMinusSSqr = kPlusDivTMinusS.multiply(k)
                .divide(TMinusS),
                TMinusSSqr = TMinusS.pow(2);
        Dfp pJ, qJ;
        SparseFieldMatrix<Dfp> tempA =
                new SparseFieldMatrix<Dfp>(koerper, k + 2, k + 2);
        /**
         * Befüllt die erste und letzte Zeile der Matrix $A$ mit den 
         * Randbedingungen.
         */
        SparseFieldVector<Dfp> tempvektor =
                new SparseFieldVector<Dfp>(koerper, k+2);
        tempvektor.setEntry(0, koerper.getOne());
        tempA.setRowVector(0, tempvektor);
        tempvektor.setEntry(0,  koerper.getZero());
        tempvektor.setEntry(k+1, koerper.getOne());
        tempA.setRowVector(k+1, tempvektor);
        /**
         * Befüllt die $k$ Zeilen der Matrix $A$, welche aus den
         * Kollokationsbedingungen hervorgehen.
         */
        for (int j = 1; j <= k; j++) {
            /**
             * Berechnung der Funktionswerte $p(x)$ und $q(x)$ für
             * $x = tau_j$.
             */
            pJ = p.value(getTau(j));
            qJ = q.value(getTau(j));
            /**
             * Berechnung des Summanden bezüglich $b_0$ der $j$-ten 
             * Kollokationsbedingung.
             */
            tempA.setEntry(j, 0, getMu(j, k - 1, true)
                    .multiply(kPlusKDivTMinusSSqr.add(getMu(j, 1, true)
                            .multiply(qJ.multiply(getMu(j, 1, true))
                            .subtract(pJ.multiply(kPlusDivTMinusS))))));
            /**
             * Berechnung des Summanden bezüglich $b_1$ der $j$-ten 
             * Kollokationsbedingung.
             */
            tempA.setEntry(j, 1, getMu(j, k-2, true).multiply(k+1)
                    .multiply((getMu(j, 1, false)
                            .multiply(k-1).subtract(getMu(j, 1, true)
                                    .multiply(koerper.getTwo()))
                            .divide(TMinusSSqr).multiply(k))
                            .add(pJ.divide(TMinusS).multiply(getMu(j, 1,false)
                                    .multiply(k+1).negate().add(1))
                                    .multiply(getMu(j, 1, true)))
                            .add(qJ.multiply(getMu(j, 2, true))
                                    .multiply(getMu(j, 1, false)))));
            /**
             * Berechnung des Summanden bezüglich $b_i, i = 2, ..., k-1$
             * der $j$-ten Kollokationsbedingung.
             */
            for (int i = 2; i < k ; i++) {
                tempA.setEntry(j, i,
                        kPlusKDivTMinusSSqr.multiply((getMu(j, 2, true)
                                .multiply(KMinusEins.getUeber(i-2)))
                                .subtract(getMu(j, 1, true)
                                        .multiply(getMu(j, 1, false))
                                        .multiply(2 * KMinusEins
                                                .getUeber(i-1)))
                                .add(getMu(j, 2, false)
                                        .multiply(KMinusEins.getUeber(i))))
                        .multiply(getMu(j, k - 1 - i, true)
                                .multiply(getMu(j, i - 2, false))).add(pJ
                                        .multiply(kPlusDivTMinusS)
                                        .multiply(getMu(j, 1, false)
                                                .multiply(KPlusEins
                                                        .getUeber(i))
                                                .negate().add(K
                                                        .getUeber(i-1)))
                                        .multiply(getMu(j, k-i, true))
                                        .multiply(getMu(j, i-1, false)))
                        .add(qJ.multiply(KPlusEins.getUeber(i)).multiply(
                                getMu(j, k+1-i, true))
                                .multiply(getMu(j, i, false))));
            }
            /**
             * Berechnung des Summanden bezüglich $b_k$ der $j$-ten 
             * Kollokationsbedingung.
             */
            tempA.setEntry(j,k, getMu(j, k-2, false).multiply(k+1).multiply(((
                    getMu(j, 1, true).multiply(k-1).subtract(getMu(j, 1,false)
                            .multiply(koerper.getTwo()))).multiply(k)
                    .divide(TMinusSSqr))
                    .add(pJ.divide(TMinusS).multiply(getMu(j, 1, false)
                            .multiply(k+1).negate().add(k))
                            .multiply(getMu(j, 1, false)))
                    .add(qJ.multiply(getMu(j, 2, false)
                            .multiply(getMu(j, 1, true)))))
                    );
            /**
             * Berechnung des Summanden bezüglich $b_(k+1)$ der $j$-ten 
             * Kollokationsbedingung.
             */
            tempA.setEntry(j, k+1, getMu(j, k - 1, false)
                    .multiply(kPlusKDivTMinusSSqr
                            .add(getMu(j, 1, false).multiply(qJ
                                    .multiply(getMu(j, 1, false))
                                    .subtract(pJ
                                            .multiply(kPlusDivTMinusS))))));
        }
//        /** Gibt die Matrix $A$ zu Testzwecken aus. */
//        for (int j = 0; j < k+2; j++) {
//            for (int i = 0; i < k+2; i++) {
//                System.out.println("A[" + j + "," + i + "] = " + 
//                tempA.getEntry(j, i));
//            }
//        }
        return tempA;
    }

    /**
     * Berechnet die Einträge der rechten Seite $v$ des zu lösenden
     * Gleichungssystems.
     * @return $v$.
     */
    private SparseFieldVector<Dfp> initialisiereV() {
        SparseFieldVector<Dfp> tempV =new SparseFieldVector<Dfp>(koerper,k+2);
        tempV.setEntry(0, eta1);
        for (int i = 1; i <= k; i++) {
            tempV.setEntry(i, f.value(getTau(i)));
        }
        tempV.setEntry(k+1, eta2);
        return tempV;
    }
    
    /**
     * Gibt $tau_j, j = 1, ..., k$ zurück.
     * @param j für das $tau_j, j = 1, ..., k$ zurückgegeben werden soll.
     * @return $tau_j$
     */
    public Dfp getTau (int j) {
        return tau[j-1];
    }
    
    /**
     * Gibt den linken Rand $s$ des Intervalls $[s, t]$ zurück.
     * @return s
     */
    private Dfp getS () {
        return xi[0];
    }
    
    /**
     * Gibt den rechten Rand $t$ des Intervalls $[s, t]$ zurück.
     * @return 
     */
    private Dfp getT () {
        return xi[xi.length-1];
    }
    
    /**
     * Gibt $mu_j^i$ oder $(1 - mu_j)^i$ zurück.
     * @param j aus $1, ..., k$ zur Auswahl von $mu_j$ oder
     * $(1 - mu_j)$.
     * @param exponent $i \in -1, 0, ..., k+1$, für den $mu_j^i$ oder
     * $(1 - mu_j)^i$.
     * @param invers {@code true}, falls $(1 - mu_j)^i$ und {@code false},
     * falls $mu_j^i$ zurückgegeben werden soll 
     * @return $mu_j^i$ oder $(1 - mu_j)^i$
     * @throws ArrayIndexOutOfBoundsException falls $i < -1$ oder
     * $i > k+1$.
     */
    public Dfp getMu (int j, int exponent, boolean invers)
            throws ArrayIndexOutOfBoundsException {
        if (exponent < -1 || exponent > k+1)
            throw new ArrayIndexOutOfBoundsException();
        if (exponent == -1)
            if (invers)
                return mu[k - j][0].reciprocal();
            else
                return mu[j - 1][0].reciprocal();
        else if (exponent == 0)
            return koerper.getOne();
        else if (invers)
                return mu[k - j][exponent-1];
            else
                return mu[j - 1][exponent-1];
        }
    
    /**
     * Gibt die bestimmte Näherungslösung zurück.
     * @return g
     */
    public BezierDarstellung getG () {
        return g;
    }

    /**
     * Gibt eine Kopie des Feldes {@code Dfp[] xi} der Gitterknoten
     * im Intervall $[s, t]$ zurück.
     * @return eine Kopie von {@code Dfp[] xi}
     */
    public Dfp[] getXi() {
        Dfp out[] = new Dfp[xi.length];
        System.arraycopy(xi, 0, out, 0, xi.length);
        return out;
    }

    /** Auszug aus den möglichen Farbnamen in Maple. */
    enum Linienfarbe {
        Brown, Crimson, Chocolate, Orange, SkyBlue, Magenta,Gold};

    /**
     * Erzeugt für ein übergebenes $k$ eine Instanz des
     * Kollokationsverfahrens zu Beispiel3 aus {G. Müllenheim, 1986}.
     * @param k die Anzahl der Kollokationspunkte.
     * @param mode Gibt an, welche Testdaten auf der Konsole ausgegeben werden
     * sollen: <p>
     * {@code mode = 1}: Gibt die Differenzen der Ableitungswerte an einigen
     * Stellen des Intervalls $[0, 1]$ aus.
     * </p><p>
     * {@code mode = 2}: Gibt die Maplebefehle zum Plotten der Näherungslösung 
     * aus. </p><p>
     * {@code mode = 3}: Gibt die Abweichungen von den Kollokations- und
     * Randbedingungen aus.
     * </p><p>
     * {@code mode = 4}: Gibt die $mu_j^i$ für alle zulässigen
     * Parameterwerte aus. </p><p>
     * {@code mode = 5}: Gibt alle $tau_j$ aus.
     * </p><p>
     * {@code mode = 6}: Gibt alle Abweichungen der Ableitungen $i$-ter
     * Ordnung, $i = 0, 1, 2$, der korrespondierenden exakten Lösung und der
     * Näherungslösung aus.
     * </p><p>
     * {@code mode = 7}: Gibt alle $xi_j$ aus.
     * </p>
     * @return Ausgabe gemäß {@code mode}.
     */
    public static String teste3 (int k, int l, int mode) {
        String ausgabe = "k = " + k + ", l = " + l + "\n";
        /** Der Körper auf dem die Funktionen definiert sind. */
        final DfpField koerper = new DfpField(genauigkeit);
        /** Der linke Rand $s$ des Kollokationsintervalls $[s, t]$. */
        final Dfp s = koerper.getOne().negate();
        /** Der rechte Rand $t$ des Kollokationsintervalls $[s, t]$. */
        final Dfp t = koerper.getOne();
        /** Der Randwert $y(s) = \eta_1$. */
        final Dfp eta1 = koerper.getZero();
        /** Der Randwert $y(t) = \eta_2$. */
        final Dfp eta2 = koerper.getZero();
        /** Die Koeffizientenfunktion $p$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getZero()}, koerper.getZero());
        /** Die Koeffizientenfunktion $q$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction q = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getTwo().negate()}, koerper.getZero());
        /** Die Funktion $f$ in $y'' + p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> f =
                new RealFieldUnivariateFunction<Dfp>() {

            public Dfp value(Dfp x) {
                Dfp tempX = x.pow(2);
                return tempX.multiply(4).multiply(tempX.exp());
            }
        };
        final BezierKollokation bkol = new BezierKollokation
                (k, l, s, t, eta1, eta2, p, q, f);
        final BezierDarstellung g = bkol.getG();
        final Beispiel3 u;
        Dfp x;
        final int n = 50;
        switch (mode) {
        case 1:
            u = new Beispiel3(koerper);
            double xD;
            final double sD = s.toDouble();
            for (int i = 0; i <= n; i++) {
                xD = sD + i * (t.toDouble() - sD)/n;
                x = koerper.newDfp(xD);
                for (int j = 0; j <= 2; j++) {
                    ausgabe += "u^(" + j + ")(" + x + ") - g^(" + j + ")(" +
                            x + ") = " + (u.getAbleitung(xD, j) -
                                    g.derivative(x, j).toDouble());
                    if (j < 2) ausgabe += "\n";
                }
                if (i < n) ausgabe += "\n";
            }
            break;
        case 2:
            String werte = "", stellen = "";
            final Dfp temp = t.subtract(s).divide(n);
            for (int i = 0; i <= n; i++) {
                x = s.add(temp.multiply(i));
                werte += g.value(x);
                stellen += x;
                if (i != n) {
                    werte += ",";
                    stellen += ",";
                }
            }
            ausgabe = "g" + k + ":=pointplot([" + stellen + " ], ["+ werte + 
                    "], legend = g^" + k + ", color = \"" + 
                    Linienfarbe.values()[(k-1)%Linienfarbe.values().length]
                            + "\", connect = true):";
            break;
        case 3:
            ausgabe += "g(" + s + ") = " + g.value(s) + "\n";
            for (int j = 1; j <= k; j++) {
                ausgabe += "tau" + j + " = " + bkol.getTau(j) + 
                        ": g'' + p * g' + q * g - f = " + 
                        (g.derivative(bkol.getTau(j), 2).add 
                                (p.value(bkol.getTau(j)).multiply 
                                (g.derivative(bkol.getTau(j), 1))).add 
                                (q.value(bkol.getTau(j)).multiply
                                (g.derivative(bkol.getTau(j), 0))).subtract 
                                (f.value(bkol.getTau(j))))
                        + "\n";
            }
            ausgabe += "g(" + t + ")= " + g.value(t);
            break;
        case 4:
            for (int j = 1; j <= k; j++) {
                for (int i = -1; i <= k+1; i++) {
                    ausgabe += "mu_" + j + "^" + i + " = " + 
                bkol.getMu(j, i, false) + "\n";
                }
            }
            break;
        case 5:
            for (int j = 1; j <= k; j++) {
                ausgabe += "tau_" + j + " = " + bkol.getTau(j) + "\n";
            }
            break;
        case 6:
            Dfp max, tempMax;
            u = new Beispiel3(koerper);
            for (int i : new int[] {0, 2}) {
                max = koerper.getZero();
                for (int j = 1; j <= k; j++) {
                    System.out.print("j = " + j + (j == k ? "\n" : ", "));
                    tempMax = g.derivative(bkol.getTau(j), i).subtract(
                            u.getAbleitung(bkol.getTau(j), i)).abs();
                    max = tempMax.greaterThan(max) ? tempMax : max;
                }
                ausgabe += "E_" + k + "^" + i + " = " + max + "\n";
            }
            break;
        case 7:
            Dfp[] tempXi = bkol.getXi();
            for (int j = 0; j <= l; j++) {
                ausgabe += "xi_" + (j + 1) + " = " + tempXi[j] + "\n";
            }
            break;
        }
        return ausgabe;
    }

    /**
     * Erzeugt für ein übergebenes $k$ eine Instanz des
     * Kollokationsverfahrens zu Beispiel4 aus {G. Müllenheim, 1986}.
     * @param k die Anzahl der Kollokationspunkte.
     * @param mode Gibt an, welche Testdaten auf der Konsole ausgegeben werden
     * sollen: <p>
     * {@code mode = 1}: Gibt die Differenzen der Ableitungswerte an einigen
     * Stellen des Intervalls $[0, 1]$ aus.
     * </p><p>
     * {@code mode = 2}: Gibt die Maplebefehle zum Plotten der Näherungslösung 
     * aus.</p><p>
     * {@code mode = 3}: Gibt die Abweichungen von den Kollokations- und
     * Randbedingungen aus.
     * </p><p>
     * {@code mode = 4}: Gibt die $mu_j^i$ für alle zulässigen
     * Parameterwerte aus. </p><p>
     * {@code mode = 5}: Gibt alle $tau_j$ aus.
     * </p><p>
     * {@code mode = 6}: Gibt alle Abweichungen der Ableitungen $i$-ter
     * Ordnung, $i = 0, 1, 2$, der korrespondierenden exakten Lösung und der
     * Näherungslösung aus.
     * </p><p>
     * {@code mode = 7}: Gibt alle $xi_j$ aus.
     * </p>
     * @return Ausgabe gemäß {@code mode}.
     */
    public static String teste4 (int k, int l, int mode) {
        String ausgabe = "k = " + k + ", l = " + l + "\n";
        /** Der Körper auf dem die Funktionen definiert sind. */
        final DfpField koerper = new DfpField(genauigkeit);
        /** Der linke Rand $s$ des Kollokationsintervalls $[s, t]$. */
        final Dfp s = koerper.getZero();
        /** Der rechte Rand $t$ des Kollokationsintervalls $[s, t]$. */
        final Dfp t = koerper.getOne();
        /** Der Randwert $y(s) = \eta_1$. */
        final Dfp eta1 = koerper.getZero();
        /** Der Randwert $y(t) = \eta_2$. */
        final Dfp eta2 = koerper.getZero();
        /** Die Koeffizientenfunktion $p$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getZero()}, koerper.getZero());
        /** Die Koeffizientenfunktion $q$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction q = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getTwo().multiply(koerper.getTwo())
                        .negate()}, koerper.getZero());
        /** Die Funktion $f$ in $y'' + p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> f = new
                RealFieldUnivariateFunction<Dfp>() {

            private final Dfp value = (koerper.getE().add(koerper.getE()
                    .reciprocal())).multiply(koerper.getTwo());
            
            public Dfp value(Dfp x) {
                return value;
            }
        };
        final BezierKollokation bkol = new BezierKollokation
                (k, l, s, t, eta1, eta2, p, q, f);
        final BezierDarstellung g = bkol.getG();
        final Beispiel4 u;
        Dfp x;
        final int n = 50;
        switch (mode) {
        case 1:
            u = new Beispiel4(koerper);
            double xD;
            final double sD = s.toDouble();
            for (int i = 0; i <= n; i++) {
                xD = sD + i * (t.toDouble() - sD)/n;
                x = koerper.newDfp(xD);
                for (int j = 0; j <= 2; j++) {
                    ausgabe += "u^(" + j + ")(" + x + ") - g^(" + j + ")(" +
                            x + ") = " + (u.getAbleitung(xD, j) -
                                    g.derivative(x, j).toDouble());
                    if (j < 2) ausgabe += "\n";
                }
                if (i < n) ausgabe += "\n";
            }
            break;
        case 2:
            String werte = "", stellen = "";
            for (int i = 0; i <= n; i++) {
                x = s.add((t.subtract(s)).divide(n).multiply(i));
                werte += g.value(x);
                stellen += x;
                if (i != n) {
                    werte += ",";
                    stellen += ",";
                }
            }
            ausgabe = "u" + k + ":=pointplot([" + stellen + " ], ["+ werte + 
                    "], legend = u^" + k + ", color = \"" + 
                    Linienfarbe.values()[(k-1)%Linienfarbe.values().length]
                            + "\", connect = true):";
            break;
        case 3:
            ausgabe += "g(" + s + ") = " + g.value(s) + "\n";
            for (int j = 1; j <= k; j++) {
                ausgabe += "tau" + j + " = " + bkol.getTau(j) + 
                        ": g'' + p * g' + q * g - f = " + 
                        (g.derivative(bkol.getTau(j), 2).add 
                                (p.value(bkol.getTau(j)).multiply 
                                (g.derivative(bkol.getTau(j), 1))).add 
                                (q.value(bkol.getTau(j)).multiply
                                (g.derivative(bkol.getTau(j), 0))).subtract 
                                (f.value(bkol.getTau(j))))
                        + "\n";
            }
            ausgabe += "g(" + t + ")= " + g.value(t);
            break;
        case 4:
            for (int j = 1; j <= k; j++) {
                for (int i = -1; i <= k+1; i++) {
                    ausgabe += "mu_" + j + "^" + i + " = " + 
                bkol.getMu(j, i, false) + "\n";
                }
            }
            break;
        case 5:
            for (int j = 1; j <= k; j++) {
                ausgabe += "tau_" + j + " = " + bkol.getTau(j) + "\n";
            }
            break;
        case 6:
            Dfp max, temp;
            u = new Beispiel4(koerper);
            for (int i : new int[] {0, 2}) {
                max = koerper.getZero();
                for (int j = 1; j <= k; j++) {
                    System.out.print("j = " + j + (j == k ? "\n" : ", "));
                    temp = g.derivative(bkol.getTau(j), i).subtract(
                            u.getAbleitung(bkol.getTau(j), i)).abs();
                    max = temp.greaterThan(max) ? temp : max;
                }
                ausgabe += "E_" + k + "^" + i + " = " + max + "\n";
            }
            break;
        case 7:
            Dfp[] tempXi = bkol.getXi();
            for (int j = 0; j <= l; j++) {
                ausgabe += "xi_" + (j + 1) + " = " + tempXi[j] + "\n";
            }
            break;
        }
        return ausgabe;
    }

    public static void main (String[] args) {        
        /**
         * Testen zweier Beispiele in verschiedenen Modi durch Ausgabe
         * interessierender Daten auf der Konsole.
         */
        for (int k = 14; k <= 14; k++) {
            for (int l = 1; l <= 1; l++) {
                System.out.println(teste3(k, l, 6));
                System.out.println(teste4(k, l, 6));
            }
        }
    }
}
