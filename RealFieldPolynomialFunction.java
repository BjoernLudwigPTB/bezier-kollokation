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
 * auf dem Datentyp $\verb!Dfp!$ angelehnt an 
 * $\verb!org.apache.commons.math3.analysis.polynomials.PolynomialFunction!$.
 */ 
public class RealFieldPolynomialFunction
        implements RealFieldUnivariateFunction<Dfp>, UnivariateFunction {

    /**
     * Das neutrale Element der Addition des Körpers, auf dem die Funktion
     * definiert ist.
     */
    private final Dfp addIdent;
    /**
     * Die Koeffizienten des Polynoms $p(x) = p_n x^n + \hdots + p_1 x + p_0,$
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
         * Der Konstruktor kopiert das übergebene Array und weist der Variable
         * $\verb!koeffizienten!$ die Kopie wie folgt zu:
         * $\mathrm{koeffizienten}[i] = p_i, i = 0, \hdots,n$.
         * @param koeffizienten Koeffizienten des Polynoms.
         * @throws NullArgumentException falls $\verb!koeffizienten null!$ ist.
         * @throws NoDataException falls $\verb!koeffizienten!$ leer ist.
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
         * Berechnet den Funktionswert des Polynoms.
         * @param x Stelle, an der die Funktion ausgewertet werden soll.
         * @return Funktionswert an der Stelle $x$.
         */
        public double value(double x) {
            return this.value(addIdent.newInstance(x)).toDouble();
        }  

        /**
         * Erzeugt aus einer $\verb!RealFieldPolynomialFunction!$ eine
         * $\verb!PolynomialFunction!$ mit entsprechend geringerer Genauigkeit.
         * @return Ein Objekt vom Typ $\verb!PolynomialFunction!$.
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
         * Ermittelt, ob die Funktion konstant ist oder nicht.
         * return $\verb!true!$, falls alle Koeffizienten $0$ sind.
         */
        public boolean istKonstantNull() {
            for (Dfp dfp : koeffizienten) {
                if (!dfp.isZero()) return false;
            }
            return true;
        } 

        /**
         * Returns the coefficients of the derivative of the polynomial with
         * the given coefficients.
         * @param coefficients Coefficients of the polynomial to differentiate.
         * @return the coefficients of the derivative or $\verb!null!$ if
         * coefficients has length $1$.
         * @throws NoDataException if $\verb!coefficients!$ is empty.
         * @throws NullArgumentException if $\verb!coefficients!$ is
         * $\verb!null!$.
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
         * Returns the derivative as a $\verb!RealFieldPolynomialFunction!$.
         * @return the derivative polynomial.
         */
        public RealFieldPolynomialFunction polynomialDerivative() {
            return new RealFieldPolynomialFunction(
                    differentiate(koeffizienten));
        }

        /**
         * Returns the derivative as a $\verb!UnivariateFunction!$.
         * @return the derivative function.
         */
        public UnivariateFunction derivative() {
            return polynomialDerivative();
        }

        /**
         * Berechnet eine Näherung der unteren Schranke des Funktionswerts
         * von $f$ im Intervall $(0, 1)$.
         * @return Eine Näherung der unteren Schranke von $f$.
         */
        public Dfp min (Dfp s, Dfp t) {
            UnivariateObjectiveFunction minF = 
                    new UnivariateObjectiveFunction(this);
            double schranke = 1e-15;
            MaxEval maxAuswertungen = new MaxEval(70);
            BrentOptimizer einfachOptimierer = new BrentOptimizer(schranke, 
                    schranke);
            MultiStartUnivariateOptimizer optimierer = new 
                    MultiStartUnivariateOptimizer(einfachOptimierer, 3, 
                            new Well19937a());
            return s.newInstance(optimierer.optimize(maxAuswertungen, 
                    minF, GoalType.MINIMIZE, new SearchInterval(s.toDouble(),
                            t.toDouble())).getValue());
        }      

        /**
         * Berechnet eine Näherung der unteren Schranke des Funktionswerts
         * von $f$ im Intervall $(0, 1)$.
         * @return Eine Näherung der unteren Schranke von $f$.
         */
        public static Dfp min (RealFieldUnivariateFunction<Dfp> f, Dfp s, Dfp t) {
            UnivariateObjectiveFunction minF = 
                    new UnivariateObjectiveFunction((UnivariateFunction) f);
            double schranke = 1e-15;
            MaxEval maxAuswertungen = new MaxEval(70);
            BrentOptimizer einfachOptimierer = new BrentOptimizer(schranke, 
                    schranke);
            MultiStartUnivariateOptimizer optimierer = new 
                    MultiStartUnivariateOptimizer(einfachOptimierer, 3, 
                            new Well19937a());
            return s.newInstance(optimierer.optimize(maxAuswertungen, 
                    minF, GoalType.MINIMIZE, new SearchInterval(s.toDouble(),
                            t.toDouble())).getValue());
        }      
}
