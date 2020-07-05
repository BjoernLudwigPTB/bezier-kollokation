import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/**
 * Enthält das Beispiel3 aus {G. Müllenheim 1986}
 */
public class Beispiel3 implements RealFieldUnivariateFunction<Dfp>,
        UnivariateFunction {
    /**
     * Der Körper auf dem die Funktion definiert ist.
     */
    private final DfpField koerper;

    /**
     * Erzeugt eine Instanz der Funktion.
     *
     * @param koerper auf dem die Funktion definiert ist.
     */
    public Beispiel3(DfpField koerper) {
        this.koerper = koerper;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double value(double x) {
        return value(koerper.newDfp(x)).toDouble();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Dfp value(Dfp x) {
        return x.pow(2).exp().subtract(
                koerper.getE().divide(koerper.getSqr2().exp()
                        .add(koerper.getSqr2().negate().exp())).multiply
                        (koerper.getSqr2().multiply(x).exp().add(koerper.getSqr2()
                                .negate().multiply(x).exp())));
    }

    /**
     * Berechnet den Wert der $\nu$-ten Ableitung der Funktion.
     *
     * @param x  Stelle, an der ausgewertet werden soll.
     * @param nu der Grad der Ableitung.
     * @return $u^{(\nu)}(x)$.
     */
    public Dfp getAbleitung(Dfp x, int nu) {
        switch (nu) {
            case 0:
                return value(x);
            case 2:
                Dfp xSqr = x.pow(2), ZweieXSqr = xSqr.exp().multiply(2);
                return ZweieXSqr.add
                        (ZweieXSqr.multiply(xSqr).multiply(2)).subtract(
                        koerper.getE().multiply(2).divide(koerper.getSqr2().exp()
                                .add(koerper.getSqr2().negate().exp())).multiply
                                (koerper.getSqr2().multiply(x).exp().add(koerper.getSqr2()
                                        .negate().multiply(x).exp())));
            default:
                return null;
        }
    }

    /**
     * Berechnet die $\nu$te Ableitung der Funktion als {@code double}.
     *
     * @param x  die Stelle an der der Ableitungswert benötigt wird.
     * @param nu der Grad der Ableitung.
     * @return $u^{(\nu)}(x)$.
     */
    public double getAbleitung(double x, int nu) {
        return getAbleitung(koerper.newDfp(x), nu).toDouble();
    }

    /**
     * Testprozedur für die Korrektheit der Implementierung
     *
     * @param args beliebige Argumente
     */
    public static void main(String[] args) {
        final DfpField koerper = new DfpField(100);
        final Dfp min = koerper.getOne().negate(), max = koerper.getOne();
        final Beispiel3 u = new Beispiel3(koerper);
        final int n = 1000;
        Dfp maxFehler = koerper.getZero(), tempFehler;
        for (int i = 0; i <= n; i++) {
            Dfp x = min.add(((max.subtract(min)).divide(n)).multiply(i));
            tempFehler = u.getAbleitung(x, 2).subtract(u.value(x)
                    .multiply(2)).subtract(x.pow(2).multiply(4)
                    .multiply(x.pow(2).exp()));
            maxFehler = tempFehler.greaterThan(maxFehler) ? tempFehler :
                    maxFehler;
        }
        System.out.println("Maximaler Fehler: " + maxFehler);
    }
}
