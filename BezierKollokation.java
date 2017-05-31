import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.*;

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
    private static final int genauigkeit = 45; 
    /**
     * Der verwendete Körper mit der Anzahl der Stellen, mit denen gerechnet
     * wird.
     */
    private DfpField koerper = new DfpField(genauigkeit);
    /** Die Anzahl der Kollokationspunkte je Teilintervall */
    private final int k;
    /** Die Kollokationspunkte $\tau_1, ..., \tau_{l \cdot k}$. */
    private final Dfp[] tau;
    /** Die Gitterknoten $s = \xi_0, \xi_1, ..., \xi_{l-1}, \xi_l = t$. */
    private final Dfp[] xi;
    /**
     * Die häufig auftretenden Faktoren $\mu_j^r := ((\tau_j-s)/(t - s))^r,$
     * $j = 1, ..., k, \quad r = 1, ..., k+1$. Es ist $mu[j-1][r-1]=\mu_j^r$.
     */
    private final Dfp[][] mu;
    /** Die Näherungslösung $g$ von $y'' + p y' + q y = f$. */
    private BezierSplineFunction g;    
    
    /**
     * Erzeugt eine Instanz des Näherungsverfahrens und stößt die Berechnung
     * der Näherungslösung an.
     * @param k Anzahl der zu verwendenden Kollokationspunkte.
     * @param l die Anzahl der Teilintervalle.
     * @param s linkes Intervallende.
     * @param t rechtes Intervallende.
     * @param eta1 linker Randwert.
     * @param eta2 rechter Randwert.
     * @param p Koeffizientenfunktion in $y'' + p y' + q y = f$.
     * @param q Koeffizientenfunktion in $y'' + p y' + q y = f$.
     * @param f in $y'' + p y' + q y = f$.
     */
    public BezierKollokation(int k, int l, Dfp s, Dfp t, Dfp eta1, Dfp eta2, 
            RealFieldUnivariateFunction<Dfp> p,
            RealFieldUnivariateFunction<Dfp> q,
            RealFieldUnivariateFunction<Dfp> f) {
        this.k = k;
        xi = initialisiereXi(l, s, t);
        tau = initialisiereTau(l);
        mu = initialisiereMus(l);
        Dfp[][] a = initialisiereA(l, p, q);
        Dfp[] v = initialisiereV(l, eta1, eta2, f);
        Dfp[] loesungBlock = new FieldBlockDecomposition(a, k, l)
                .solve(v);
        BezierFunction[] functions = new BezierFunction[l];
        for (int i = 0; i < l; i++) {
            functions[i] = new BezierFunction(koerper, ArrayUtils
                    .subarray(loesungBlock, i * (k + 2), (i + 1) * (k + 2)),
                    getXi(i), getXi(i+1));            
        }
        g = new BezierSplineFunction(getXi(), functions);
    }
    
    /**
     * Berechnet alle $\mu_j^r := \mu(\tau_j)^r$ und
     * $(1 - \mu_j)^r := (1 - \mu(\tau_j))^r$ und speichert diese. Vor
     * dem Aufruf dieser Prozedur muss das Feld {@code Dfp[] tau} mit den
     * Kollokationspunkten initialisiert sein. Die Berechnung der
     * Funktionswerte ist dabei auf ein Minimum beschränkt unter Ausnutzung
     * von <p>
     * $\mu_j = 1 - \mu_{k-j-1}$,
     * </p>
     * und äquivalent <p>
     * $1 - \mu_j = \mu_{k-j-1}$.</p>
     * @param l die Anzahl der Teilintervalle.
     * @return ein Feld {@code Dfp[] mu} mit mu[j][r] $= \mu_j^r$.
     */ 
    private Dfp[][] initialisiereMus (int l) {
        Dfp[][] tempMu = new Dfp[l*k][k+1];
        for (int i = 0; i < l; i++) {
            BezierFunction tempG = new BezierFunction(koerper, 
                    new Dfp[] {koerper.getZero(), koerper.getZero()},
                    getXi(i), getXi(i+1));
            int index1 = i * k;
            for (int j = 0; j < k; j++) {
                int index2 = index1 + j;
                if (j > k - j)
                    tempMu[index2][0] = koerper.getOne()
                    .subtract(tempMu[index1 + k - j - 1][0]);
                else
                    tempMu[index2][0] = tempG.getMu(getTau(index2 + 1));
                for (int r = 1; r <= k; r++) {
                    tempMu[index2][r] = tempMu[index2][r-1]
                            .multiply(tempMu[index2][0]);
                }
            }
        }
        return tempMu;
    }
    
    /**
     * Berechnet die streng monoton steigende Folge
     * $\tau_j, j = 1, ..., l \cdot k$ im Intervall $[s, t]$.
     */
    private Dfp[] initialisiereTau (int l) {
        Dfp[] temptau = new Dfp[l*k];
        GaussLegendrePunkte rhos = new GaussLegendrePunkte(k, koerper);
        for (int i = 0; i < l; i++) {
            Dfp tempPlus = getXi(i).add(getXi(i+1)), tempMinus = getXi(i+1)
                    .subtract(getXi(i));
            for (int j = 0; j < k; j++) {
                temptau[(i * k) + j] = tempPlus.add(tempMinus.multiply(rhos
                        .getRho(j))).divide(koerper.getTwo());
            }
        }
        return temptau;
    }
    
    /**
     * Berechnet eine streng monoton steigende Folge
     * $\xi_i, i = 0, ..., l$ im Intervall $[s, t]$ mit $\xi_0 = s,\xi_l = t$.
     * @param s das linke Intervallende von $[s, t]$.
     * @param t das rechte Intervallende von $[s, t]$.
     * @param l die Anzahl der Teilintervalle.
     */
    private Dfp[] initialisiereXi (int l, Dfp s, Dfp t) {
        Dfp[] tempXi = new Gitterknoten(l, s, t).getXi();
        return tempXi;
    }
    
    /**
     * Berechnet die Einträge der Koeffizientenmatrix $A$ des zu lösenden
     * Gleichungssystems.
     * @param l die Anzahl der Teilintervalle.
     * @param p in $y'' + p y' + q y = f$.
     * @param q in $y'' + p y' + q y = f$.
     * @return $A$.
     */
    private Dfp[][] initialisiereA (int l,
            RealFieldUnivariateFunction<Dfp> p,
            RealFieldUnivariateFunction<Dfp> q) {
        /*
         * Für $k = 1$ entsteht bei der blockweisen Erstellung die
         * Besonderheit, dass die Blöcke der Länge $k + 2 = 3$ nicht mehr 
         * die Stetigkeitsbedingungen aufnehmen können, da sich diese
         * jeweils auf vier Bézierpunkte beziehen. Deshalb wird in diesem
         * Fall die Matrix geringfügig modifiziert.
         */
        Dfp[][] tempA = new Dfp[l * (k + 2)][(k == 1 && l > 1) ? k + 3 : k+2];
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
            tempA[tempA.length-1][i] = koerper.getZero();
        }
        tempA[tempA.length-1][tempA[0].length - 1] = koerper.getOne();
        /* Erzeugen der Blockstruktur durch verschachtelte Schleifen. */
        for (int i = 0; i < l; i++) {
            {
                /*
                 * Variablen zum Zwischenspeichern wiederholt benötigter
                 * Zwischenergebnisse.
                 */
                Dfp deltaXiMinus = deltaXi;
                deltaXi = getXi(i+1).subtract(getXi(i));
                /*
                 * Befüllt die $2$ Zeilen der Matrix $A$, welche aus den
                 * Stetigkeitsbedingungen für den Übergang vom $i-1$-ten zum
                 * $i$-ten Teilintervall hervorgehen, $i = 1, \hdots, l-1$.
                 */
                if (i > 0) {
                    /* Die Stetigkeitsbedingung bezüglich der $C^1$
                     * -Stetigkeit. */
                    tempA[i*(k+2)-1][0] = deltaXi;
                    tempA[i*(k+2)-1][1] = deltaXiMinus.add(deltaXi).negate();
                    tempA[i*(k+2)-1][2] = koerper.getZero();
                    tempA[i*(k+2)-1][3] = deltaXiMinus;
                    /* Die Stetigkeitsbedingung bezüglich der $C^0$
                     * -Stetigkeit. */
                    tempA[i*(k+2)][0] = koerper.getZero();
                    tempA[i*(k+2)][1] = koerper.getOne();
                    tempA[i*(k+2)][2] = koerper.getOne().negate();
                    tempA[i*(k+2)][3] = koerper.getZero();
                    for (int j = 4; j < tempA[0].length; j++) {
                        tempA[i*(k+2)-1][j] = koerper.getZero();
                        tempA[i*(k+2)][j] = koerper.getZero();
                    }
                }
            }
            Dfp kPlusDivDeltaXi = deltaXi.reciprocal().multiply(k + 1),
                    kPlusTimesKDivDeltaXiSqr = kPlusDivDeltaXi.multiply(k)
                    .divide(deltaXi), deltaXiSqr = deltaXi.pow(2);
            /*
             * Befüllt die $k$ Zeilen der Matrix $A$, welche aus den
             * Kollokationsbedingungen des $i$-ten Teilintervalls hervorgehen.
             */
            for (int j = 1; j <= k; j++) {
                /*
                 * Berechnung der Funktionswerte $p(\tau_{ik+j})$ und
                 * $q(\tau_{ik+j})$ für den $j$-ten Kollokationspunkt des 
                 * $i$-ten Teilintervalls.
                 */
                Dfp pJ = p.value(getTau(i, j)), qJ = q.value(getTau(i, j));
                /*
                 * Berechnung des Summanden bezüglich $b_{i0}$ der $j$-ten 
                 * Kollokationsbedingung.
                 */
                tempA[i*(k+2)+j][0] = getMu(i, j, k - 1, true)
                        .multiply(kPlusTimesKDivDeltaXiSqr.add(getMu(i, j, 1,
                                true).multiply(qJ.multiply(getMu(i, j, 1,
                                        true)).subtract(pJ.multiply(
                                                kPlusDivDeltaXi)))));
                /*
                 * Berechnung des Summanden bezüglich $b_{i1}$ der $j$-ten 
                 * Kollokationsbedingung.
                 */
                tempA[i*(k+2)+j][1] = getMu(i, j, k-2, true)
                        .multiply(k+1).multiply((getMu(i, j, 1, false)
                                .multiply(k-1).subtract(getMu(i, j, 1, true)
                                        .multiply(koerper.getTwo()))
                                .divide(deltaXiSqr).multiply(k))
                                .add(pJ.divide(deltaXi).multiply(getMu(i, j, 
                                        1,false).multiply(k+1).negate().add(
                                                koerper.getOne()))
                                        .multiply(getMu(i, j, 1, true)))
                                .add(qJ.multiply(getMu(i, j, 2, true))
                                        .multiply(getMu(i, j, 1, false))));
                /*
                 * Berechnung des Summanden bezüglich
                 * $b_{i\kappa}, \kappa = 2, ..., k-1$ der $j$-ten
                 * Kollokationsbedingung.
                 */
                for (int kappa = 2; kappa < k ; kappa++) {
                    tempA[i * (k + 2) + j][kappa] =
                            kPlusTimesKDivDeltaXiSqr.multiply((getMu(i, j, 2,
                                    true).multiply(binomM.getUeber(kappa-2)))
                                    .subtract(getMu(i, j, 1, true)
                                            .multiply(getMu(i, j, 1, false))
                                            .multiply(2 * 
                                                    binomM.getUeber(kappa-1)))
                                    .add(getMu(i, j, 2, false).multiply(binomM
                                            .getUeber(kappa))))
                            .multiply(getMu(i, j, k - 1 - kappa, true)
                                    .multiply(getMu(i, j, kappa - 2, false)))
                            .add(pJ.multiply(kPlusDivDeltaXi)
                                    .multiply(getMu(i, j, 1, false)
                                            .multiply(binomP.getUeber(kappa))
                                            .negate()
                                            .add(binomK.getUeber(kappa-1)))
                                    .multiply(getMu(i, j, k-kappa, true))
                                    .multiply(getMu(i, j, kappa-1, false)))
                            .add(qJ.multiply(binomP.getUeber(kappa))
                                    .multiply(getMu(i, j, k+1-kappa, true))
                                    .multiply(getMu(i, j, kappa, false)));
                }
                /*
                 * Berechnung des Summanden bezüglich $b_{ik}$ der $j$-ten 
                 * Kollokationsbedingung.
                 */
                tempA[i*(k+2)+j][k] = getMu(i, j, k-2, false)
                        .multiply(k+1).multiply(((getMu(i, j, 1, true)
                                .multiply(k-1)
                                .subtract(getMu(i, j, 1,false)
                                        .multiply(koerper.getTwo())))
                                .multiply(k).divide(deltaXiSqr)).add(pJ
                                        .divide(deltaXi)
                                        .multiply(getMu(i, j, 1, false)
                                                .multiply(k+1).negate()
                                                .add(k))
                                        .multiply(getMu(i, j, 1, false)))
                                .add(qJ.multiply(getMu(i, j, 2, false)
                                        .multiply(getMu(i, j, 1, true)))));
                /*
                 * Berechnung des Summanden bezüglich $b_{i,k+1}$ der $j$-ten 
                 * Kollokationsbedingung.
                 */
                tempA[i*(k+2)+j][k+1] = getMu(i,j,k-1, false)
                        .multiply(kPlusTimesKDivDeltaXiSqr
                                .add(getMu(i, j, 1, false).multiply(qJ
                                        .multiply(getMu(i, j, 1, false))
                                        .subtract(pJ.multiply(
                                                kPlusDivDeltaXi)))));
                if (k == 1 && l > 1) tempA[i*(k+2)+j][k+2] = koerper.getZero();
            }
        }
//        /* Gibt die Matrix $A$ zu Testzwecken aus. */
//        for (int i = 0; i < l * (k+2); i++) {
//            for (int j = 0; j < l * (k+2); j++) {
//                System.out.println("A[" + i + "," + j + "] = " + 
//                tempA.getEntry(i, j));
//            }
//        }
        return tempA;
    }

    /**
     * Berechnet die Einträge der rechten Seite $v$ des zu lösenden
     * Gleichungssystems.
     * @param l die Anzahl der Teilintervalle.
     * @param eta1 Randwert $y(s) = \eta_1$.
     * @param eta2 Randwert $y(t) = \eta_2$.
     * @param f in $y'' + p y' + q y = f$.
     * @return $v$.
     */
    private Dfp[] initialisiereV (int l, Dfp eta1,
            Dfp eta2, RealFieldUnivariateFunction<Dfp> f) {
        Dfp[] tempV = new Dfp[l * (k + 2)];
        tempV[0] = eta1;
        for (int i = 0; i < l; i++) {
            int index = i * (k+2);
            if (i > 0) {
                tempV[index - 1] = koerper.getZero();
                tempV[index] = koerper.getZero();
            }
            for (int j = 1; j <= k; j++) {
                tempV[index + j] = f.value(getTau(i, j));
            }
        }
        tempV[l * (k+2) - 1] = eta2;
//        /* Gibt den Vektor $v$ zu Testzwecken aus. */
//        for (int j = 0; j < l * (k+2); j++) {
//                System.out.println("v[" + j + "] = " + 
//                        tempV[j]);
//        }
        return tempV;
    }
    
    /**
     * Gibt $\tau_j, j = 1, ..., l \cdot k$ zurück.
     * @param j für das $\tau_j, j = 1, ..., l \cdot k$ zurückgegeben werden
     * soll.
     * @return $\tau_j$
     */
    public Dfp getTau (int j) {
        return tau[j-1];
    }
    
    /**
     * Gibt $\tau_j, j = 1, ..., k$ aus dem $i$-ten Teilintervall zurück.
     * @param i für das $i$-te Teilintervall.
     * @param j für das $\tau_j, j = 1, ..., k$ zurückgegeben werden soll.
     * @return $\tau_{ik+j}$
     */
    public Dfp getTau (int i, int j) {
        return getTau(i*k + j);
    }
    
    /**
     * alte Version von getMu, die nur die neue mit einem Standardwert
     * $i = 0$ für das $0$-te Teilintervall aufruft. */
    public Dfp getMu (int j, int exponent, boolean invers)
            throws ArrayIndexOutOfBoundsException {
        return getMu(0, j, exponent, invers);
    }
    
    /**
     * Gibt $\mu_j^{\operatorname{exponent}}$ oder
     * $(1 - \mu_j)^{\operatorname{exponent}}$ zurück.
     * @param i das $i$-te Teilintervall $i = 0, \hdots, l - 1$.
     * @param j aus $1, ..., k$ zur Auswahl von $\mu_j$ oder
     * $(1 - \mu_j)$.
     * @param exponent ${\operatorname{exponent}} \in \{0, ..., k+1\}$, für
     * den $\mu_j^{\operatorname{exponent}}$ oder
     * $(1 - \mu_j)^{\operatorname{exponent}}$.
     * @param invers {@code true}, falls
     * $(1 - \mu_j)^{\operatorname{exponent}}$ und {@code false},
     * falls $\mu_j^{\operatorname{exponent}}$ zurückgegeben werden soll. 
     * @return $\mu_j^{\operatorname{exponent}}$ oder
     * $(1 - \mu_j)^{\operatorname{exponent}}$
     * @throws ArrayIndexOutOfBoundsException falls 
     * ${\operatorname{exponent}} < 0$ oder ${\operatorname{exponent}} > k+1$.
     */
    public Dfp getMu (int i, int j, int exponent, boolean invers)
            throws ArrayIndexOutOfBoundsException {
        if (exponent < -1 || exponent > k + 1 || i < 0 
                || i > xi.length + 1)
            throw new ArrayIndexOutOfBoundsException();
        if (exponent == -1)
            if (invers)
                return mu[(i+1) * k - j][0].reciprocal();
            else
                return mu[i * k + j - 1][0].reciprocal();
        else if (exponent == 0)
            return koerper.getOne();
        else if (invers) {
            return mu[(i+1) * k - j][exponent-1];
        }
        else {
            return mu[i * k + j - 1][exponent-1];
        }
    }
    
    /**
     * Gibt die bestimmte Näherungslösung zurück.
     * @return g
     */
    public BezierSplineFunction getG () {
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
    
    /**
     * Gibt $\xi_i, i = 0, ..., l$ zurück.
     * @param i für das $\xi_i, i = 1, ..., l$ zurückgegeben werden soll.
     * @return $\xi_i$
     */
    public Dfp getXi (int i) {
        return xi[i];
    }

    /** Auszug aus den möglichen Farbnamen in Maple. */
    enum Linienfarbe {
        Brown, Crimson, Chocolate, Orange, SkyBlue, Magenta,Gold};

    /**
     * Erzeugt für ein übergebenes $k$ und ein $l$ eine Instanz des
     * Kollokationsverfahrens zu Beispiel3 aus {G. Müllenheim, 1986}.
     * @param k die Anzahl der Kollokationspunkte.
     * @param l Anzahl der Teilintervalle.
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
     * {@code mode = 4}: Gibt die $\mu_j^i$ für alle zulässigen
     * Parameterwerte aus. </p><p>
     * {@code mode = 5}: Gibt alle $\tau_j$ aus.
     * </p><p>
     * {@code mode = 6}: Gibt alle Abweichungen der Ableitungen $i$-ter
     * Ordnung, $i = 0, 1, 2$, der korrespondierenden exakten Lösung und der
     * Näherungslösung aus.
     * </p><p>
     * {@code mode = 7}: Gibt alle $\xi_j$ aus.
     * </p>
     * @return Ausgabe gemäß {@code mode}.
     */
    public static String teste3 (int k, int l, int mode) {
        String ausgabe = "k = " + k + ", l = " + l + "\n";
        /* Der Körper auf dem die Funktionen definiert sind. */
        final DfpField koerper = new DfpField(genauigkeit);
        /* Der linke Rand $s$ des Kollokationsintervalls $[s, t]$. */
        final Dfp s = koerper.getOne().negate();
        /* Der rechte Rand $t$ des Kollokationsintervalls $[s, t]$. */
        final Dfp t = koerper.getOne();
        /* Der Randwert $y(s) = \eta_1$. */
        final Dfp eta1 = koerper.getZero();
        /* Der Randwert $y(t) = \eta_2$. */
        final Dfp eta2 = koerper.getZero();
        /* Die Koeffizientenfunktion $p$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getZero()});
        /* Die Koeffizientenfunktion $q$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction q = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getTwo().negate()});
        /* Die Funktion $f$ in $y'' + p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> f =
                new RealFieldUnivariateFunction<Dfp>() {

            public Dfp value(Dfp x) {
                Dfp tempX = x.pow(2);
                return tempX.multiply(4).multiply(tempX.exp());
            }
        };
        final BezierKollokation bkol = new BezierKollokation
                (k, l, s, t, eta1, eta2, p, q, f);
        final BezierSplineFunction g = bkol.getG();
        final Beispiel3 u;
        Dfp x;
        final int n = 200;
        switch (mode) {
        case 1:
            u = new Beispiel3(koerper);
            double xD;
            final double sD = s.toDouble();
            for (int i = 0; i <= n; i++) {
                xD = sD + i * (t.toDouble() - sD)/n;
                x = koerper.newDfp(xD);
                for (int j : new int[] {0, 2}) {
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
            ausgabe = "g" + k + l + ":=pointplot([" + stellen + " ], ["+ werte + 
                    "], legend = \"g" + k + l + "\", color = \"" + 
                    Linienfarbe.values()[(l+k-1)%Linienfarbe.values().length]
                            + "\", connect = true):";
            break;
        case 3:
            ausgabe += "g(" + s + ") = " + g.value(s) + "\n";
            for (int j = 1; j <= l * k; j++) {
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
            for (int j = 1; j <= l * k; j++) {
                for (int i = 0; i <= k+1; i++) {
                    ausgabe += "mu_" + j + "^" + i + " = " + 
                bkol.getMu(j, i, false) + "\n";
                }
            }
            break;
        case 5:
            for (int i = 0; i < l; i++) {
                for (int j = 1; j <= k; j++) {
                    ausgabe += "tau_" + (i * k + j) + " = " + bkol.getTau(i,
                            j) + "\n";
                }
            }
            break;
        case 6:
            u = new Beispiel3(koerper);
            Dfp max = koerper.getZero();
            for (int j = 1; j <= l * k; j++) {
//                Ausgabe auf der Konsole während der Berechnungen, um ein 
//                Bild von der Dauer zu bekommen.
//                System.out.print("j = " + j + (j == k * l ? "\n" : ", "));
                Dfp temp = g.value(bkol.getTau(j)).subtract(
                        u.value(bkol.getTau(j))).abs();
                if (temp.greaterThan(max)) max = temp;
            }
            ausgabe += "E_\\Delta^" + k + (Integer.toString(k).length() == 1 ?
                    " " : "") + " = " + max + "\n";
            max = koerper.getZero();
            for (int i = 0; i <= l; i++) {
//                Ausgabe auf der Konsole während der Berechnungen, um ein
//                Bild von der Dauer zu bekommen.
//                System.out.print("i = " + i + (i == l ? "\n" : ", "));
               Dfp temp = g.value(bkol.getXi(i)).subtract(
                        u.value(bkol.getXi(i))).abs();
                if (temp.greaterThan(max)) max = temp;
            }
            ausgabe += "E_\\xi^" + l;
            for (int m = -4; m < Integer.toString(k).length() - Integer
                    .toString(l).length() - (Integer.toString(k)
                            .length() == 1 ? 0 : 1); m++) {
                ausgabe += " ";
                }
            ausgabe += " = " + max + "\n";
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
     * Erzeugt für ein übergebenes $k$ und ein $l$ eine Instanz des
     * Kollokationsverfahrens zu Beispiel4 aus {G. Müllenheim, 1986}.
     * @param k die Anzahl der Kollokationspunkte.
     * @param l Anzahl der Teilintervalle.
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
     * {@code mode = 4}: Gibt die $\mu_j^i$ für alle zulässigen
     * Parameterwerte aus. </p><p>
     * {@code mode = 5}: Gibt alle $\tau_j$ aus.
     * </p><p>
     * {@code mode = 6}: Gibt alle Abweichungen der Ableitungen $i$-ter
     * Ordnung, $i = 0, 1, 2$, der korrespondierenden exakten Lösung und der
     * Näherungslösung aus.
     * </p><p>
     * {@code mode = 7}: Gibt alle $\xi_j$ aus.
     * </p>
     * @return Ausgabe gemäß {@code mode}.
     */
    public static String teste4 (int k, int l, int mode) {
        String ausgabe = "k = " + k + ", l = " + l + "\n";
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
        /* Die Koeffizientenfunktion $p$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getZero()});
        /* Die Koeffizientenfunktion $q$ in $y'' + p y' + q y = f$. */
        final RealFieldPolynomialFunction q = new RealFieldPolynomialFunction(
                new Dfp[] {koerper.getTwo().multiply(koerper.getTwo())
                        .negate()});
        /* Die Funktion $f$ in $y'' + p y' + q y = f$. */
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
        final BezierSplineFunction g = bkol.getG();
        final Beispiel4 u;
        Dfp x;
        final int n = 200;
        switch (mode) {
        case 1:
            u = new Beispiel4(koerper);
            double xD;
            final double sD = s.toDouble();
            for (int i = 0; i <= n; i++) {
                xD = sD + i * (t.toDouble() - sD)/n;
                x = koerper.newDfp(xD);
                for (int j : new int[] {0, 2}) {
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
            ausgabe = "u" + k + l + ":=pointplot([" + stellen + " ], ["+ werte + 
                    "], legend = \"u" + k + l + "\", color = \"" + 
                    Linienfarbe.values()[(l+k-1)%Linienfarbe.values().length]
                            + "\", connect = true):";
            break;
        case 3:
            ausgabe += "g(" + s + ") = " + g.value(s) + "\n";
            for (int j = 1; j <= l * k; j++) {
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
            for (int j = 1; j <= l * k; j++) {
                for (int i = 0; i <= k+1; i++) {
                    ausgabe += "mu_" + j + "^" + i + " = " + 
                bkol.getMu(j, i, false) + "\n";
                }
            }
            break;
        case 5:
            for (int i = 0; i < l; i++) {
                for (int j = 1; j <= k; j++) {
                    ausgabe += "tau_" + (i * k + j) + " = " + bkol.getTau(i,
                            j) + "\n";
                }
            }
            break;
        case 6:
            u = new Beispiel4(koerper);
            Dfp max = koerper.getZero();
            for (int j = 1; j <= l * k; j++) {
//                Ausgabe auf der Konsole während der Berechnungen, um ein 
//                Bild von der Dauer zu bekommen.
//                System.out.print("j = " + j + (j == k * l ? "\n" : ", "));
                Dfp temp = g.value(bkol.getTau(j)).subtract(
                        u.value(bkol.getTau(j))).abs();
                if (temp.greaterThan(max)) max = temp;
            }
            ausgabe += "E_\\Delta^" + k + (Integer.toString(k).length() == 1 ?
                    " " : "") + " = " + max + "\n";
            max = koerper.getZero();
            for (int i = 0; i <= l; i++) {
//                Ausgabe auf der Konsole während der Berechnungen, um ein
//                Bild von der Dauer zu bekommen.
//                System.out.print("i = " + i + (i == l ? "\n" : ", "));
                Dfp temp = g.value(bkol.getXi(i)).subtract(
                        u.value(bkol.getXi(i))).abs();
                if (temp.greaterThan(max)) max = temp;
            }
            ausgabe += "E_\\xi^" + l;
            for (int m = -4; m < Integer.toString(k).length() - Integer
                    .toString(l).length() - (Integer.toString(k)
                            .length() == 1 ? 0 : 1); m++) {
                ausgabe += " ";
                }
            ausgabe += " = " + max + "\n";
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
     * Testen zweier Beispiele in verschiedenen Modi durch Ausgabe
     * interessierender Daten auf der Konsole.
     */
    public static void main (String[] args) {
        for (int l : new int[] {2, 3, 4}) {
            for (int k : new int[] {1, 3, 5, 10}) {
//                System.out.println("Beispiel 3:");
                System.out.println(teste3(k, l, 6));
//                System.out.println("Beispiel 4:");
//                System.out.println(teste4(k, l, 6));
            }
        }
    }
}
