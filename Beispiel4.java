import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/**
 * Enthält das Beispiel4 aus {G. Müllenheim 1986}
 */
public class Beispiel4 implements RealFieldUnivariateFunction<Dfp>,
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
    public Beispiel4(DfpField koerper) {
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
        return (koerper.getTwo().multiply(x).subtract(koerper.getOne()))
                .cosh().subtract(koerper.getOne().cosh());
    }

    /**
     * Berechnet den Wert der $\nu$-ten Ableitung der Funktion.
     *
     * @param x  Stelle, an der ausgewertet werden soll.
     * @param nu der Grad der Ableitung.
     * @return $u^{(\nu)}(x)$.
     */
    public Dfp getAbleitung(Dfp x, int nu) {
        return switch (nu) {
            case 0 -> value(x);
            case 1 -> (x.multiply(2).subtract(1)).sinh().multiply(2);
            case 2 -> (x.multiply(2).subtract(1)).cosh().multiply(4);
            default -> null;
        };
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
        final Dfp min = koerper.getZero();
        final Dfp max = koerper.getOne();
        final int n = 10;
        final Beispiel4 u = new Beispiel4(koerper);
        for (int i = 0; i <= n; i++) {
            Dfp x = min.add(max.subtract(min).multiply(i)).divide(n);
            System.out.println("u''(" + x + ") - 4 * u(" + x +
                    ") - 4 * cosh(1) = " + (u.getAbleitung(x, 2)
                    .subtract(u.value(x).multiply(4))
                    .subtract((koerper.getOne()).cosh()
                            .multiply(4))));
        }
    }
}
