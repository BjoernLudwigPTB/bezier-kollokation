import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

/**
 * Implementierung eines Polynoms in Bézierdarstellung angelehnt an 
 * {@code org.apache.commons.math3.analysis.polynomials.PolynomialFunction}.
 */
public class BezierDarstellung implements RealFieldUnivariateFunction<Dfp> {

     /**
      *  Der Körper auf dem die Bézierdarstellung implementiert wird.
      */
        DfpField koerper;
        /**
         * The bezierpoints of the polynomial, ordered by degree -- i.e.,
         * $b[0]$ is $b_0$ and $b[n]$ is the coefficient of
         * ${}_s^t B_n^n$ where n is the degree of the polynomial.
         */
        private final Dfp[] b;
        /**
         * Die Intervallgrenzen $[s, t]$ über dem die Basispolynome 
         * definiert sind.
         */
        private final Dfp s, t;
        
        /**
         * Construct a polynomial with the given bezierpoints. The first
         * element of {@code Dfp[] b} is the bezierpoint for $B_0^{(k+1)}$. 
         * Higher degree bezierpoints follow in sequence. The degree of the 
         * resulting polynomial is the index of the last element of the array.
         * <p>
         * The constructor makes a copy of the input array and assigns the
         * copy to {@code Dfp[] b}.
         * </p>
         * @param b bezierpoints.
         * @throws NullArgumentException if {@code b} is {@code null}.
         * @throws NoDataException if {@code b} is empty.
         */
        public BezierDarstellung(DfpField koerper, Dfp[] b, Dfp s, Dfp t)
                throws NullArgumentException, NoDataException {
            this.koerper = koerper;
            this.s = s;
            this.t= t;
            MathUtils.checkNotNull(b);
            int n = b.length;
            if (n == 0) {
                throw new NoDataException(
                        LocalizedFormats
                        .EMPTY_POLYNOMIALS_COEFFICIENTS_ARRAY);
            }
            this.b = new Dfp[n];
            System.arraycopy(b, 0, this.b, 0, n);
        }
        
        /** 
         * Auswertung des Polynoms in Bézierdarstellung an der Stelle $x$.
         * @param x die Stelle an der ausgewertet werden soll.
         * @return Funktionswert an der Stelle $x$.
         */
        public Dfp value (Dfp x) {
            return deCasteljau(x, mu(x, s, t), b.length - 1, 0);
        }
        
        /** 
         * Auswertung der $m$-ten Ableitung, $m = 0, 1, 2$ des
         * Polynoms in Bézierdarstellung an der Stelle $x$.
         * @param x die Stelle an der ausgewertet werden soll.
         * @param m Grad der Ableitung, die ausgewertet werden soll.
         * @return Ableitungsswert an der Stelle $x$.
         */
        public Dfp derivative (Dfp x, int m) {
            Dfp tempValue = koerper.getZero();
            Dfp mu = mu(x, s, t);
            switch (m) {
            case 0:
                tempValue = value(x);
                break;
            case 1:
                tempValue = (t.subtract(s)).reciprocal()
                .multiply(b.length - 1).multiply(deCasteljau(x, mu, b.length
                        -2, 1).subtract(deCasteljau(x, mu, b.length - 2, 0)));
                break;
            case 2:
                tempValue = (t.subtract(s).pow(2)).reciprocal()
                .multiply(b.length - 1).multiply(b.length - 2)
                .multiply(deCasteljau(x, mu, b.length - 3, 0)
                        .subtract(koerper.getTwo().multiply(
                                deCasteljau(x, mu, b.length - 3, 1)))
                        .add(deCasteljau(x, mu, b.length - 3, 2)));
                break;
            }
            return tempValue;
        }

        /**
         * Returns a copy of the bezierpoints.
         * <p>
         * Changes made to the returned copy will not affect the bezierpoints
         * of the polynomial.
         * </p>
         * @return a fresh copy of the bezierpoints array.
         */
        public Dfp[] getBezierpunkte() {
            return (Dfp[]) b.clone();
        }
        
        /**
         * Berechnet den häufig auftretenden Faktor $\mu(x) := (x-s)/(t-s)$.
         * @param x Wert, für den $\mu(x)$ berechnet werden soll.
         * @return $\mu(x)$
         */
        public static Dfp mu (Dfp x, Dfp s, Dfp t) {
            return (x.subtract(s)).divide(t.subtract(s));
        }
        
        /**
         * Führt den Algorithmus von deCasteljau durch.
         * @param x Stelle die ausgewertet werden soll.
         * @param mu Quotient $(x-s)/(t-s)$ für die Berechnung der Rekursion
         * mit den Intervallgrenzen $[s, t]$.
         * @param r oberer Index $b_i^r$.
         * @param i unterer Index $b_i^r$.
         * @param b Feld der Bézier-Punkte.
         * @return Für $i = 0, r = n}$ der Funktionswert des Polynoms in
         * Bézier-Darstellung mit Bézier-Punkten {@code doube[] b} an der
         * Stelle $x$.
         */
        public Dfp deCasteljau (Dfp x, Dfp mu, int r, int i) {
            if (r == 0)
                return b[i];
            else
                return  (mu.multiply(deCasteljau(x, mu, r - 1, i + 1)))
                        .add(koerper.getOne().subtract(mu)
                                .multiply(deCasteljau(x, mu, r - 1, i)));
        }
    }