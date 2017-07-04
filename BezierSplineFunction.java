import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.exception.*;
import org.apache.commons.math3.exception.util.LocalizedFormats;

/**
 * Represents a Bézier spline function, most of the implementation being taken
 * from
 * $\verb!apache.commons.math3.analysis.polynomials.PolynomialSplineFunction!$.
 * A $\textbf{Bézier spline function}$ consists of a set of
 * $\textit{BezierFunction}$s and an ascending array of domain
 * $\textit{gridknots}$, determining the intervals over which the spline
 * function is defined by the constituent polynomials.
 * The domain of the polynomial spline function is
 * $[s := \text{smallest knot}, \text{largest knot} =: t]$. Attempts to
 * evaluate the spline at values smaller than $s$ will give the function value
 * of the first polynomial and analoguely for values bigger than $t$ the
 * function value of the last polynomial.
 * The value of the polynomial spline function for an argument $x$ is computed
 * as follows:
 * 1. The knot array is searched to find the segment to which $x$ belongs.
 * If $x$ is less than the smallest knot point the value returned is
 * $\operatorname{functions}[0](x)$.
 * 2. Let $j$ be the index of the largest knot point that is less
 * than or equal to $x$. The value returned is
 * $\operatorname{functions}[j](x)$.
 */
public class BezierSplineFunction implements RealFieldUnivariateFunction<Dfp>{
    /** Spline segment interval delimiters of size $l + 1$ for $l$ segments.*/
    private final Dfp[] knoten;
    /**
     * The polynomial functions that make up the spline.  The first element
     * determines the value of the spline over the first subinterval, the
     * second over the second, etc..
     */
    private final BezierFunction[] functions;
    /**
     * Number of spline segments. It is equal to the number of polynomials and
     * to the number of partition points $- 1$.
     */
    private final int l;

    /**
     * Construct a polynomial spline function with the given segment
     * delimiters and polynomials. The constructor copies both arrays and
     * assigns the copies to the $\verb!gitterKnoten!$ and
     * $\verb!polynomials!$ properties, respectively.
     * @param knoten Spline segment interval delimiters.
     * @param functions Polynomial functions that make up the spline.
     * @throws NullArgumentException if either of the input arrays is
     * $\verb!null!$.
     * @throws NumberIsTooSmallException if knoten has length less than 2.
     * @throws DimensionMismatchException if 
     * $\verb!functions.length \!= knoten.length - 1!$.
     *
     */
    public BezierSplineFunction (Dfp knoten[], BezierFunction functions[])
        throws NullArgumentException, NumberIsTooSmallException,
               DimensionMismatchException{
        if (knoten == null ||
            functions == null) {
            throw new NullArgumentException();
        }
        if (knoten.length < 2) {
            throw new NumberIsTooSmallException(
                    LocalizedFormats.NOT_ENOUGH_POINTS_IN_SPLINE_PARTITION,
                    2, knoten.length, false);
        }
        if (knoten.length - 1 != functions.length) {
            throw new DimensionMismatchException(functions.length,
                    knoten.length);
        }
        this.l = knoten.length -1;
        this.knoten = new Dfp[l + 1];
        System.arraycopy(knoten, 0, this.knoten, 0, l + 1);
        this.functions = new BezierFunction[l];
        System.arraycopy(functions, 0, this.functions, 0, l);
    }

    /**
     * Compute the value for the function.
     * See above for details on the algorithm for
     * computing the value of the function.
     * @param x Point for which the function value should be computed.
     * @return the value.
     */
    @Override
    public Dfp value (Dfp x) {
        return functions[getInterval(x)].value(x);
    }

    /**
     * Compute the derivative for the function.
     * See above for details on the algorithm for
     * computing the derivative of the function.
     * @param x Point for which the function derivative should be computed.
     * @return the value.
     */
    public Dfp derivative (Dfp x, int nu) {
        return functions[getInterval(x)].derivative(x, nu);
    }

    /**
     * Get the number of spline segments It is also the number of functions
     * and the number of knoten $- 1$.
     * @return the number of spline segments.
     */
    public int getL() {
        return l;
    }

    /**
     * Get a copy of the Bézierfunctions array.
     * It returns a fresh copy of the array. Changes made to the copy will
     * not affect the $\verb!functions!$ property.
     * @return the functions.
     */
    public BezierFunction[] getPolynomials() {
        return functions.clone();
    }

    /**
     * Get an array copy of the Gitterknoten.
     * It returns a fresh copy of the array. Changes made to the copy
     * will not affect the knots property.
     * @return the Gitterknoten.
     */
    public Dfp[] getKnots() {
        return knoten.clone();
    }
    
    /**
     * Berechnet das zugehörige Intervall für $x$.
     * @param x der Wert für den das Intervall bestimmt werden soll.
     * @return $i \in \{0, \hdots, l - 1\}$.
     */
    private int getInterval (Dfp x) {
        if (x.lessThan(knoten[1]) || knoten.length == 2) {
            return 0;
        } else {
            int i;
            for (i = 1; i < l - 1; i++) {
                if (x.lessThan(knoten[i+1])) break;
            }
            return i;
        }
    }
}
