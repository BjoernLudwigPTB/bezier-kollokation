import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

/**
 * Repräsentiert eine Polynomfunktion beliebigen Grades basierend
 * auf dem Datentyp {@code Dfp} angelehnt an 
 * {@code org.apache.commons.math3.analysis.polynomials.PolynomialFunction}.
 */ 
public class RealFieldPolynomialFunction
        implements RealFieldUnivariateFunction<Dfp>, UnivariateFunction {

    /**
     * Das neutrale Element der Addition des Körpers, auf dem die Funktion
     * definiert ist.
     */
    private final Dfp addIdent;
    /**
     * Die Koeffizienten des Polynoms
     * <p>
     * $p(x) = p_n x^n + \hdots + p_1 x + p_0,$
     * $\mathrm{koeffizienten}[i] = p_i, \quad i = 0, \hdots ,n.$
     * </p>
     */
    private final Dfp[] koeffizienten;

        /**
         * Konstruiert ein Polynom $p(x) = p_n x^n + \hdots + p_1 x + p_0$
         * mit den übergebenen Koeffizienten. Das erste Element des Arrays
         * ist der konstante Term. Höhergradige Koeffizienten folgen in Reihe.
         * Der Grad des sich ergebenden Polynoms ist der Index des letzten von
         * Null verschiedenen Elements des Array, falls nicht alle Elemente
         * Null sind.
         * <p>
         * Der Konstruktor kopiert das übergebene Array und weist der Variable
         * {@code koeffizienten} die Kopie wie folgt zu:</p>
         * <p>
         * $\mathrm{koeffizienten}[i] = p_i, i = 0, \hdots,n$</p>
         *
         * @param koeffizienten Koeffizienten des Polynoms.
         * @throws NullArgumentException falls {@code koeffizienten}
         * {@code null} ist.
         * @throws NoDataException falls {@code koeffizienten} leer ist.
         */
        public RealFieldPolynomialFunction(Dfp[] koeffizienten, Dfp addIdent)
                throws NullArgumentException, NoDataException {
            super();
            this.addIdent = addIdent;
            MathUtils.checkNotNull(koeffizienten);
            int j = koeffizienten.length;
            if (j == 0) {
                throw new NoDataException(LocalizedFormats
                        .EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
            }
            while ((j > 1) && (koeffizienten[j - 1] == addIdent)) {
                --j;
            }
            this.koeffizienten = new Dfp[j];
            System.arraycopy(koeffizienten, 0, this.koeffizienten, 0, j);
        }

        /**
         * @param x Stelle, an der die Funktion ausgewertet werden soll.
         * @return Funktionswert an der Stelle $x$.
         */
        public Dfp value (Dfp x) {
            Dfp wert = addIdent;
            for (int i = 0; i < koeffizienten.length; i++) {
                wert = wert.add((x.pow(i)).multiply(koeffizienten[i]));
            } 
            return wert;
        }

        /**
         * Berechnet den Funktionswert des Polynoms als {@code double}.
         * @param x Stelle, an der die Funktion ausgewertet werden soll.
         * @return Funktionswert an der Stelle $x$.
         */
        public double value(double x) {
            DfpField koerper = addIdent.getField();
            return koerper.newDfp(this.value(koerper.newDfp(x))).toDouble();
        }  

        /**
         * Erzeugt aus einer {@code RealFieldPolynomialFunction} eine
         * {@code PolynomialFunction} mit entsprechend geringerer Genauigkeit.
         * @param f ist die umzuwandelnde {@code RealFieldPolynomialFunction}.
         * @return Ein Objekt vom Typ {@code PolynomialFunction}.
         */
        public PolynomialFunction getPolynomialFunction () {
            int anzahl = koeffizienten.length;
            double[] doubleKoeffizienten = new double[anzahl];
            for (int i = 0; i < anzahl; i++) {
                doubleKoeffizienten[i] = addIdent.getField()
                        .newDfp(koeffizienten[i]).toDouble();
            }
            return new PolynomialFunction(doubleKoeffizienten);
        }
}
