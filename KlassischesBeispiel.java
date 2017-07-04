import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/** Repräsentiert ein Beispiel aus {$\textsc{de Boor, Schwartz}$ 1973}. */
public class KlassischesBeispiel implements RealFieldUnivariateFunction<Dfp>,
UnivariateFunction {
    /** Der Körper auf dem die Funktion definiert ist. */
    private final DfpField koerper;
    
    /**
     * Erzeugt eine Instanz der Funktion.
     * @param koerper auf dem die Funktion definiert ist.
     */
    public KlassischesBeispiel (DfpField koerper) {
        this.koerper = koerper;
    }
    
    /** {@inheritDoc} */
    @Override
    public double value(double x) {
        return value(koerper.newDfp(x)).toDouble();
    }
    
    /** {@inheritDoc} */
    @Override
    public Dfp value(Dfp x) {
        return (koerper.getTwo().multiply(x).subtract(koerper.getOne()))
                .cosh().subtract(koerper.getOne().cosh());
    }
    
    /**
     * Berechnet den Wert der $\nu$-ten Ableitung der Funktion.
     * @param x Stelle, an der ausgewertet werden soll.
     * @param nu der Grad der Ableitung.
     * @return $u^{(\nu)}(x)$.
     */
    public Dfp getAbleitung(Dfp x, int nu) {
        switch (nu) {
        case 0:
            return value(x);
        case 1:
            return koerper.getTwo().multiply(x).subtract(koerper.getOne())
                    .sinh().multiply(koerper.getTwo());
        case 2:
            return koerper.getTwo().multiply(x).subtract(koerper.getOne())
                    .cosh().multiply(koerper.getTwo().pow(2));
        default:
            return null;
        }
    }
    
    /** 
     * Berechnet die $\nu$te Ableitung der Funktion als $\verb!double!$.
     * @param x die Stelle an der der Ableitungswert benötigt wird.
     * @param nu der Grad der Ableitung.
     * @return $u^{(\nu)}(x)$.
     */
    public double getAbleitung(double x, int nu) {
        return getAbleitung(koerper.newDfp(x), nu).toDouble();
    }
}
