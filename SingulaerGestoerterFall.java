import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/**
 * Repräsentiert den singulär gestörten Fall und behandelt das Verfahren auf
 * den verschiedenen Gittern.
 */
public class SingulaerGestoerterFall extends BezierKollokation {

    /**
     * Erzeugt eine Instanz im singulär gestörten Fall. Dabei wird der
     * allgemeine Konstruktor unverändert aufgerufen, aber die ermittelte
     * Näherungslösung in den statischen Methoden der Klasse anders behandelt.
     * @param k Anzahl der Kollokationspunkte je Gitterintervall.
     * @param xi das Gitter auf dem die Näherungslösung berechnet werden soll.
     * @param epsilon singulärer Störungsparameter.
     * @param eta1 linker Randwert.
     * @param eta2 rechter Randwert.
     * @param p Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param q Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param f in $-\varepsilon y'' - p y' + q y = f$.
     */
    public SingulaerGestoerterFall(int k, Gitter<Dfp> xi, Dfp epsilon,
            Dfp eta1, Dfp eta2, RealFieldPolynomialFunction p,
            RealFieldUnivariateFunction<Dfp> q,
            RealFieldUnivariateFunction<Dfp> f) {
        super(k, xi, epsilon, eta1, eta2, p, q, f);
    }
    
    /**
     * Erzeugt für ein übergebenes $k$ und ein Feld $l$ Instanzen des
     * Kollokationsverfahrens für das Reaktionsdiffusionsproblem und wertet
     * diese in Abhängigkeit von $\verb!mode!$ aus.
     * @param k die Anzahl der Kollokationspunkte je Gitterintervall.
     * @param l Feld mit den interessierenden Anzahlen von Gitterintervallen.
     * @param mode Gibt an, welche Testdaten auf der Konsole ausgegeben werden
     * sollen:
     * $\verb!mode = 1!$: Berechnet die Fehler und die experimentelle
     * Konvergenzordnung für das Reaktionsdiffusionsproblem.
     * $\verb!mode = 2!$: Gibt die Maplebefehle zum Plotten der Näherungslösung aus.
     * $\verb!mode = 3!$: Gibt die Abweichungen von den Kollokations- und
     * Randbedingungen aus.
     * @return Ausgabe gemäß $\verb!mode!$ in Form eines formatierten $\verb!String!$.
     */
    public static String berechneReaktion (int k, int[] l, final Dfp epsilon,
            int mode) {
        /* Der linke Rand $s$ des Kollokationsintervalls $[s, t]$. */
        final Dfp s = epsilon.getZero();
        /* Der rechte Rand $t$ des Kollokationsintervalls $[s, t]$. */
        final Dfp t = epsilon.getOne();
        /* Der Randwert $y(s) = \eta_1$. */
        final Dfp eta1 = epsilon.getZero();
        /* Der Randwert $y(t) = \eta_2$. */
        final Dfp eta2 = epsilon.getZero();
        /* Die Koeffizientenfunktion $p$ in $-\varepsilon^2 y'' - p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction (
                new Dfp[] {epsilon.getZero()});
        /* Die Koeffizientenfunktion $q$ in $-\varepsilon^2 y'' - p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> q = new
                RealFieldUnivariateFunction<Dfp>() {
            
            public Dfp value(Dfp x) {
                return x.cos().add(x.multiply(x)).add(x.getOne());
            }
        };
        /* Die Funktion $f$ in $-\varepsilon^2 y'' - p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> f = new
                RealFieldUnivariateFunction<Dfp>() {
            public Dfp value(Dfp x) {
                return x.pow(9d/2).add(x.sin());
            }
        };
        /* Häufig benötigte Ausdrücke. */
        Dfp qu_0 = epsilon.getTwo().pow(2).reciprocal(),
                qu_1 = qu_0,
                sigma_0 = epsilon.getTwo().multiply(epsilon.getTwo()),
                sigma_1 = epsilon.getTwo().multiply(epsilon.getTwo());
        String ausgabe = "k = " + k + ", ";
        switch (mode) {
        case 1:
            /*
             * $\verb!max[i]!$ enthält die maximale Abweichung zur Vergleichslösung
             * bei $\verb!l[i]!$ Gitterintervallen.
             */
            Dfp[] max = new Dfp[l.length];
            for (int i = 0; i < l.length; i++) {
                max[i] = epsilon.getZero();
            }
            ausgabe += "epsilon = 10^" + ((int) (epsilon.log().toDouble()
                    /Math.log(10)))+ ": " + "\n";
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> originalGitter = new ShishkinGitter(l[i],s,t,qu_0,
                        qu_1, sigma_0, sigma_1, epsilon.getTwo(), epsilon);
                max[i] = berechneAbweichung(k, originalGitter,
                        epsilon.pow(2), eta1, eta2, p, q, f);
                ausgabe += "l = "; 
                for (int m = 0; m < 3 - Integer.toString(l[i]).length(); m++){
                    ausgabe += " "; 
                };
                ausgabe += l[i] + ": E_\\xi = " + max[i] + "\n";
                /*
                 * Berechnet die experimentelle Konvergenzordnung wie in der
                 * Ausarbeitung angegeben.
                 */
                if (i > 0) {
                        for (int m = 0; m < 2-Integer.toString(l[i]).length();
                                m++) {
                            ausgabe += " ";
                        }
                        ausgabe += "   \\alpha_" + l[i] + " = " +max[i].divide
                                (max[i-1]).log().divide(Math.log(Math.log(l[i]
                                        ) / ( 2*Math.log(l[i-1])))) + "\n";
                }
                ausgabe += "\n"; 
            }            
            break;
        case 2:
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> gitter = new ShishkinGitter(l[i], s, t, qu_0,
                        qu_1, sigma_0, sigma_1, epsilon.getTwo(), epsilon);
                final SingulaerGestoerterFall kollokation = new 
                        SingulaerGestoerterFall(k, gitter, epsilon.pow(2), 
                                eta1, eta2, p, q, f);
                final BezierSplineFunction g = kollokation.getG();
                String werte = "", stellen = "";
                for (int j = 0; j < l[i]; j++) {
                    werte += g.value(kollokation.getXi(j)).toDouble();
                    stellen += kollokation.getXi(j).toDouble();
                    werte += ",";
                    stellen += ",";
                }
                werte += g.value(kollokation.getXi(l[i])).toDouble();
                stellen += kollokation.getXi(l[i]).toDouble();
                ausgabe = "u" + k + "_" + l[i] + "_" + Math.abs((int) (epsilon
                        .log().toDouble() / Math.log(10))) + ":=pointplot([" +
                        stellen + " ], [" + werte + "], legend = \"u" + k + 
                        "_" + l[i] + "_" + Math.abs((int) (epsilon.log()
                                .toDouble() / Math.log(10))) +"\", color = \""
                        +Linienfarbe.values()[i%(Linienfarbe.values().length)]
                                + "\", connect = true):";
            }
            break;
        case 3:
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> gitter = new ShishkinGitter(l[i], s, t, qu_0,
                        qu_1, sigma_0, sigma_1, epsilon.getTwo(), epsilon);               
                        final SingulaerGestoerterFall kollokation = new 
                        SingulaerGestoerterFall(k, gitter, epsilon.pow(2), 
                                eta1, eta2, p, q, f);
                final BezierSplineFunction g = kollokation.getG();
                ausgabe += "g(" + s + ") = " + g.value(s) + "\n";
                for (int j = 1; j <= l[i] * k; j++) {
                    Dfp temp = kollokation.getTau(j);
                    ausgabe += "tau" + j + " = " + temp.toDouble() + 
                            ": -epsilon g'' - p * g' + q * g - f = " + 
                            (epsilon.pow(2).negate().multiply(
                                    g.derivative(temp, 2)).subtract(p
                                            .value(temp)
                                            .multiply(g.derivative(temp, 1)))
                                    .add(q.value(temp).multiply
                                            (g.derivative(temp, 0)))
                                    .subtract(f.value(temp)))
                            + "\n";
                }
                ausgabe += "g(" + t + ")= " + g.value(t);
            }
        }
        return ausgabe;
    }
    
    /**
     * Erzeugt für die übergebenen Daten eine Vergleichslösung auf einem 7-mal
     * feineren Gitter und berechnet die maximale Abweichung der
     * Funktionswerte an den Gitterknoten des Originals.
     * @param k Anzahl der Kollokationspunkte je Gitterintervall.
     * @param xi das Gitter auf dem die Originalnäherung berechnet wurde.
     * @param epsilon singulärer Störungsparameter.
     * @param eta1 linker Randwert.
     * @param eta2 rechter Randwert.
     * @param p Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param q Koeffizientenfunktion in $-\varepsilon y'' - p y' + q y = f$.
     * @param f in $-\varepsilon y'' - p y' + q y = f$.
     */
    private static Dfp berechneAbweichung (int k, Gitter<Dfp> originalGitter,
            Dfp epsilon, Dfp eta1, Dfp eta2, RealFieldPolynomialFunction p,
            RealFieldUnivariateFunction<Dfp> q,
            RealFieldUnivariateFunction<Dfp> f) {
        Dfp max = epsilon.getZero(); 
        Gitter<Dfp> modifiziertesGitter = new ModifiziertesGitter(
                originalGitter,7);
        final SingulaerGestoerterFall originalKollokation = new 
                SingulaerGestoerterFall(k, originalGitter, epsilon,
                        eta1, eta2, p, q, f), modifizierteKollokation = new
                        SingulaerGestoerterFall(k,modifiziertesGitter, 
                                epsilon, eta1, eta2, p, q, f);
        final BezierSplineFunction originalG = originalKollokation
                .getG(),modifiziertesG=modifizierteKollokation.getG();
        for (int j = 0; j <= originalGitter.getXi().length-1; j++) {
            Dfp temp = originalG.value(originalGitter.getXi(j))
                    .subtract(modifiziertesG.value(originalGitter
                            .getXi(j))).abs();
            if (temp.greaterThan(max)) {
                max = temp;
            }
        }
        return max;
    }
    
    /**
     * Erzeugt für ein übergebenes $k$ und ein Feld $l$ Instanzen des
     * Kollokationsverfahrens für das Konvektionsdiffusionsproblem und wertet
     * diese in Abhängigkeit von $\verb!mode!$ aus.
     * @param k die Anzahl der Kollokationspunkte je Gitterintervall.
     * @param l Feld mit den interessierenden Anzahlen von Gitterintervallen.
     * @param mode Gibt an, welche Testdaten auf der Konsole ausgegeben werden
     * sollen:
     * $\verb!mode = 11!$: Berechnet die Fehler und die experimentelle
     * Konvergenzordnung für das Konvektionsdiffusionsproblem auf einem
     * $\textsc{Shishkin}$-Gitter.
     * $\verb!mode = 12!$: Berechnet die Fehler und die experimentelle
     * Konvergenzordnung für das Konvektionsdiffusionsproblem auf einem
     * $\textsc{Bakhvalov}$-Gitter.
     * $\verb!mode = 2!$: Gibt die Maplebefehle zum Plotten der Näherungslösung  aus.
     * $\verb!mode = 3!$: Gibt die Abweichungen von den Kollokations- und
     * Randbedingungen aus.
     * @return Ausgabe gemäß $\verb!mode!$ in Form eines formatierten $\verb!String!$.
     */
    public static String berechneKonvektion (int k, int[] l, final
            Dfp epsilon, int mode) {
        /* Der linke Rand $s$ des Kollokationsintervalls $[s, t]$. */
        final Dfp s = epsilon.getZero();
        /* Der rechte Rand $t$ des Kollokationsintervalls $[s, t]$. */
        final Dfp t = epsilon.getOne();
        /* Der Randwert $y(s) = \eta_1$. */
        final Dfp eta1 = epsilon.getZero();
        /* Der Randwert $y(t) = \eta_2$. */
        final Dfp eta2 = epsilon.getZero();
        /* Die Koeffizientenfunktion $p$ in $-\varepsilon y'' - p y' + q y = f$. */
        final RealFieldPolynomialFunction p = new RealFieldPolynomialFunction(
                new Dfp[] {epsilon.getOne()});
        /* Die Koeffizientenfunktion $q$ in $-\varepsilon y'' - p y' + q y = f$. */
        final RealFieldPolynomialFunction q = new RealFieldPolynomialFunction(
                new Dfp[] {epsilon.getTwo()});
        /* Die Funktion $f$ in $-\varepsilon y'' - p y' + q y = f$. */
        final RealFieldUnivariateFunction<Dfp> f = new
                RealFieldUnivariateFunction<Dfp>() {
            public Dfp value(Dfp x) {
                return x.subtract(epsilon.getOne()).exp();
            }
        };
        /* Häufig benötigter Ausdruck. */
        Dfp qu = epsilon.getTwo().reciprocal(), sigma = epsilon.getOne();
        Dfp[] max;
        String ausgabe = "k = " + k + ", ";
        switch (mode) {
        case 11:
            /*
             * $\verb!max[i]!$ enthält die maximale Abweichung zur Vergleichslösung
             * bei $\verb!l[i]!$ Gitterintervallen.
             */
            max = new Dfp[l.length];
            for (int i = 0; i < l.length; i++) {
                max[i] = epsilon.getZero();
            }
            ausgabe += "epsilon = 10^" + ((int) (Math.log(epsilon.toDouble())
                    /Math.log(10)))+ ": " + "\n";
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> originalGitter = new ShishkinGitter(l[i], s, t,
                        qu, sigma, epsilon, epsilon);
                max[i] = berechneAbweichung(k, originalGitter, epsilon, 
                        eta1, eta2, p, q, f);
                ausgabe += "l = "; 
                for (int m = 0; m < 3 - Integer.toString(l[i]).length(); m++){
                    ausgabe += " "; 
                };
                ausgabe += l[i] + ": E_\\xi = " + max[i] + "\n";
                if (i > 0) {
                        for (int m = 0; m < 2-Integer.toString(l[i]).length();
                                m++) {
                            ausgabe += " ";
                        };
                        ausgabe += "    \\alpha_" + l[i] + " = " + 
                                max[i].divide(max[i-1]).log().divide(
                                        Math.log(Math.log(l[i])/(2*
                                                Math.log(l[i-1])))) +"\n";
                }
                ausgabe += "\n"; 
            }            
            break;
        case 12:
            /*
             * $\verb!max[i]!$ enthält die maximale Abweichung zur Vergleichslösung
             * bei $\verb!l[i]!$ Gitterintervallen.
             */
            max = new Dfp[l.length];
            for (int i = 0; i < l.length; i++) {
                max[i] = epsilon.getZero();
            }
            ausgabe += "epsilon = 10^" + ((int) (Math.log(epsilon.toDouble())
                    /Math.log(10)))+ ": " + "\n";
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> originalGitter = new BakhvalovGitter(l[i], s, t,
                        qu, sigma, epsilon, epsilon);
                max[i] = berechneAbweichung(k, originalGitter, epsilon, 
                        eta1, eta2, p, q, f);
                ausgabe += "l = "; 
                for (int m = 0; m < 3 - Integer.toString(l[i]).length(); m++){
                    ausgabe += " "; 
                };
                ausgabe += l[i] + ": E_\\xi = " + max[i] + "\n";
                if (i > 0) {
                        for (int m = 0; m < 2-Integer.toString(l[i]).length();
                                m++) {
                            ausgabe += " ";
                        };
                        ausgabe += "    \\alpha_" + l[i] + " = " + 
                                max[i-1].divide(max[i]).log().divide(
                                        Math.log(2)) + "\n";
                }
                ausgabe += "\n"; 
            }            
            break;
        case 2:
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> gitter = new BakhvalovGitter(l[i], s, t, qu, 
                        sigma, epsilon, epsilon);
                final SingulaerGestoerterFall kollokation = new 
                        SingulaerGestoerterFall(k, gitter, epsilon, eta1,
                                eta2, p, q, f);
                final BezierSplineFunction g = kollokation.getG();
                String werte = "", stellen = "";
                for (int j = 0; j < l[i]; j++) {
                    werte += g.value(kollokation.getXi(j)).toDouble();
                    stellen += kollokation.getXi(j).toDouble();
                    werte += ",";
                    stellen += ",";
                }
                werte += g.value(kollokation.getXi(l[i])).toDouble();
                stellen += kollokation.getXi(l[i]).toDouble();
                ausgabe = "g" + k + "_" + l[i] + "_" + Math.abs((int) (epsilon
                        .log().toDouble() / Math.log(2))) + ":=pointplot([" + 
                        stellen + " ], [" + werte + "], legend = \"u" + k + 
                        "_" + l[i] + "_" + Math.abs((int) (epsilon.log()
                                .toDouble() / Math.log(10)))+"\", color = \""+
                                Linienfarbe.values()[i]+"\", connect =true):";
            }
            break;
        case 3:
            for (int i = 0; i < l.length; i++) {
                Gitter<Dfp> gitter = new BakhvalovGitter(l[i], s, t, qu,sigma,
                        epsilon, epsilon);               
                        final SingulaerGestoerterFall kollokation = new 
                        SingulaerGestoerterFall(k, gitter, epsilon, eta1,
                                eta2, p, q, f);
                final BezierSplineFunction g = kollokation.getG();
                ausgabe += "g(" + s + ") = " + g.value(s) + "\n";
                for (int j = 1; j <= l[i] * k; j++) {
                    Dfp temp = kollokation.getTau(j);
                    ausgabe += "tau" + j + " = " + temp.toDouble() + 
                            ": -epsilon g'' - p * g' + q * g - f = " + 
                            (epsilon.negate().multiply(g.derivative(temp, 2))
                                    .subtract(p.value(temp)
                                            .multiply(g.derivative(temp, 1)))
                                    .add(q.value(temp).multiply
                                            (g.derivative(temp, 0)))
                                    .subtract(f.value(temp)))
                            + "\n";
                }
                ausgabe += "g(" + t + ")= " + g.value(t);
            }
        }
        return ausgabe;
    }
    
    /**
     * Berechnen verschiedener Daten in Abhängigkeit des angegebenen Modus
     * und Ausgabe auf der Konsole.
     */
    public static void main (String[] args) {
        DfpField koerper = new DfpField(BezierKollokation.getGenauigkeit());
        int kleinstesL = 3, größtesL = 6, deltaL = größtesL - kleinstesL;
        int[] l = new int[deltaL + 1];
        l[0] = (int) Math.pow(2, kleinstesL);
        for (int i = 1; i <= deltaL; i++) {
            l[i] = l[i-1] * 2;
        }
        Dfp[] epsilon = new Dfp[] {koerper.newDfp(10).pow(-24),
                koerper.newDfp(10).pow(-36), koerper.newDfp(10).pow(-48)};
        for (int i = 0; i < epsilon.length; i++) {
            for (int k : new int[] {1, 2, 4}) {
              System.out.println(berechneReaktion(k, l, epsilon[i], 1));
//                System.out.println(berechneKonvektion(k, l, epsilon[i],11));
            }      
        }
    }
}
