import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.MathUtils;

/**
 * Repräsentiert ein Polynom in Bézierdarstellung angelehnt an 
 * {@code org.apache.commons.math3.analysis.polynomials.PolynomialFunction}.
 */
public class BezierFunction implements RealFieldUnivariateFunction<Dfp> {

     /** Der Körper auf dem die Bézierdarstellung implementiert wird. */
        DfpField koerper;
        /**
         * The bezierpoints of the polynomial, ordered by degree -- i.e.,
         * {@code b[0]} is $b_0$ and {@code b[n]} is the coefficient of
         * ${}_s^t B_n^n$ where $n$ is the degree of all the Bernstein
         * polynomials the function is composed of.
         */
        private final Dfp[] b;
        /**
         * Die Intervallgrenzen $[s, t]$ über dem die Basispolynome 
         * definiert sind und häufig auftretende Rechenschritte.
         */
        private final Dfp tMinusS, tMinusSSqr, SDivTMinusS;
        
        /**
         * Construct a polynomial with the given bezierpoints. The first
         * element of {@code Dfp[] b} is the bezierpoint for ${}_s^t B_0^{k+1}$. 
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
        public BezierFunction(DfpField koerper, Dfp[] b, Dfp s, Dfp t)
                throws NullArgumentException, NoDataException {
            this.koerper = koerper;
            tMinusS = t.subtract(s);
            tMinusSSqr = tMinusS.pow(2);
            SDivTMinusS = s.divide(tMinusS);
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
            return deCasteljau(x, getMu(x), b.length - 1, 0);
        }
        
        /** 
         * Auswertung der $m$-ten Ableitung, $m = 0, 1, 2, \hdots$ des
         * Polynoms in Bézierdarstellung an der Stelle $x$.
         * @param x die Stelle an der ausgewertet werden soll.
         * @param nu Grad der Ableitung, die ausgewertet werden soll.
         * @return Ableitungsswert an der Stelle $x$.
         */
        public Dfp derivative (Dfp x, int nu) {
            Dfp mu = getMu(x);
            if (nu == 0) return value(x);
            /* 
             * Switch-Anweisung, weil die ersten drei Ableitungen
             * $\nu = 0, 1, 2$ aus Geschwindigkeitsgründen hart codiert sind.
             */
            switch (nu) {
            case 0:
                return value(x);
            case 1:
                return (tMinusS).reciprocal()
                .multiply(b.length - 1).multiply(deCasteljau(x, mu, b.length
                        -2, 1).subtract(deCasteljau(x, mu, b.length - 2, 0)));
            case 2:
                return (tMinusSSqr).reciprocal()
                .multiply(b.length - 1).multiply(b.length - 2)
                .multiply(deCasteljau(x, mu, b.length - 3, 0)
                        .subtract(koerper.getTwo().multiply(
                                deCasteljau(x, mu, b.length - 3, 1)))
                        .add(deCasteljau(x, mu, b.length - 3, 2)));
            default:
                Dfp tempValue = koerper.getZero();
                Binomialkoeffizient mUeberI = new Binomialkoeffizient(nu);
                int n = 1;
                for (int i = 0; i <= nu; i++) {
                    tempValue = tempValue.add(((nu-i)%2 == 0) ? 
                            deCasteljau(x, mu, b.length-1-nu, i)
                            .multiply(mUeberI.getUeber(i))
                            : deCasteljau(x, mu, b.length-1-nu, i)
                            .multiply(mUeberI.getUeber(i)).negate());
                    if (i < nu) n *= b.length - 1 - i; 
                }
                return tempValue.multiply(n).divide(tMinusS.pow(nu));
            }
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
        
        /**
         * Berechnet den häufig auftretenden Faktor $\mu(x) := \frac{x-s}{t-s}$.
         * @param x Wert, für den $\mu(x)$ berechnet werden soll.
         * @return $\mu(x)$
         */
        public Dfp getMu (Dfp x) {
            return x.divide(tMinusS).subtract(SDivTMinusS);
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
         * Gibt den Körper zurück, auf dem die Funktion definiert ist.
         * @return a fresh copy of the bezierpoints array.
         */
        public DfpField getKoerper() {
            return koerper;
        }
        
        /**
         * Berechnet den häufig auftretenden Faktor $\mu(x) := (x-s)/(t-s)$.
         * @param x Wert, für den $\mu(x)$ berechnet werden soll.
         * @return $\mu(x)$
         */
        public static Dfp mu (Dfp x, Dfp s, Dfp t) {
            return (x.subtract(s)).divide(t.subtract(s));
        }        
    }