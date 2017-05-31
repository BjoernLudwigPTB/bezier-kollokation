import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;
import org.apache.commons.math3.random.Well19937a;
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
         * @param koeffizienten Koeffizienten des Polynoms.
         * @throws NullArgumentException falls {@code koeffizienten}
         * {@code null} ist.
         * @throws NoDataException falls {@code koeffizienten} leer ist.
         */
        public RealFieldPolynomialFunction(Dfp[] koeffizienten)
                throws NullArgumentException, NoDataException {
            super();
            this.addIdent = koeffizienten[0].getZero();
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
            return this.value(addIdent.newInstance(x)).toDouble();
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
                doubleKoeffizienten[i] = addIdent
                        .newInstance(koeffizienten[i]).toDouble();
            }
            return new PolynomialFunction(doubleKoeffizienten);
        }

        /**
         * Berechnet eine Näherung der unteren Schranke des Funktionswerts
         * im Intervall $(0, 1)$.
         * @return Eine Näherung der unteren Schranke von {@code f}.
         */
        public Dfp min () {
            UnivariateFunction f = new Beispiel4(addIdent.getField());
            UnivariateObjectiveFunction minF = 
                    new UnivariateObjectiveFunction(f);
            double schranke = 1.0e-15;
            MaxEval maxAuswertungen = new MaxEval(40);
            BrentOptimizer einfachOptimierer = new BrentOptimizer(schranke, 
                    schranke);
            MultiStartUnivariateOptimizer optimierer = new 
                    MultiStartUnivariateOptimizer(einfachOptimierer, 1, 
                            new Well19937a());
            System.out.println(einfachOptimierer.getMax());
            return addIdent.newInstance(optimierer.optimize(maxAuswertungen, 
                    minF, GoalType.MINIMIZE, new SearchInterval(-1, 1))
                    .getValue() - schranke);
        }       

        /**
         * Returns the coefficients of the derivative of the polynomial with
         * the given coefficients.
         * @param coefficients Coefficients of the polynomial to differentiate.
         * @return the coefficients of the derivative or {@code null} if
         * coefficients has length $1$.
         * @throws NoDataException if {@code coefficients} is empty.
         * @throws NullArgumentException if {@code coefficients} is
         * {@code null}.
         */
        protected static Dfp[] differentiate(Dfp[] coefficients)
            throws NullArgumentException, NoDataException {
            MathUtils.checkNotNull(coefficients);
            int n = coefficients.length;
            if (n == 0) {
                throw new NoDataException(LocalizedFormats
                        .EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
            }
            if (n == 1) {
                return new Dfp[]{coefficients[0].getField().getZero()};
            }
            Dfp[] result = new Dfp[n - 1];
            for (int i = n - 1; i > 0; i--) {
                result[i - 1] = coefficients[i].multiply(i);
            }
            return result;
        }

        /**
         * Returns the derivative as a {@link RealFieldPolynomialFunction}.
         * @return the derivative polynomial.
         */
        public RealFieldPolynomialFunction polynomialDerivative() {
            return new RealFieldPolynomialFunction(
                    differentiate(koeffizienten));
        }

        /**
         * Returns the derivative as a {@link UnivariateFunction}.
         * @return the derivative function.
         */
        public UnivariateFunction derivative() {
            return polynomialDerivative();
        }
}
