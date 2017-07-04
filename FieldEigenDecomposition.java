import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.util.LocalizedFormats;

/**
 * Calculates the eigenvalues of a real matrix taking the algortihm from
 * $\verb!EigenDecomposition!$ and adapting it to $\verb!Dfp!$ matrices. This class is similar
 * in spirit to the $\verb!EigenvalueDecomposition!$ class from the JAMA library, with the
 * limitation to calculating Eigenvalues of symmetric tridiagonal matrices
 * with only zeros on the main diagonal to save operations. This
 * implementation is based on the paper by $\textsc{A. Drubrulle, R.S. Martin}$ and
 * $\textsc{J.H. Wilkinson}$ "The Implicit QL Algorithm" in $\textsc{Wilksinson}$ and $\textsc{Reinsch}$
 * (1971) Handbook for automatic computation, vol. 2, Linear algebra,
 * Springer-Verlag, New-York.
 */
public class FieldEigenDecomposition {
    /** Der Körper über dem die Zerlegung berechnet werden soll. */
    private final DfpField koerper;
    /** Maximum number of iterations accepted. */
    private final static byte maxIter = 30;
    /** The eigenvalues. */
    private final Dfp[] eigenvalues;

    /**
     * Calculates the eigen decomposition of the symmetric tridiagonal
     * matrix.
     * @param neben of the tridiagonal form.
     * @throws MaxCountExceededException if the algorithm fails to converge.
     * @since 3.1
     */ 
    public FieldEigenDecomposition(final Dfp[] neben) {
            this.koerper = neben[0].getField();
            eigenvalues = findEigenValues(neben);
        }
    
    /**
     * Find eigenvalues ($\textsc{Dubrulle}$ et al., 1971)
     * @param neben Secondary diagonal of the tridiagonal matrix.
     */
    private Dfp[] findEigenValues(Dfp[] neben) {
        final int n = neben.length;
        Dfp[] tempEigenvalues = new Dfp[n];
        final Dfp[] e = neben;
        for (int i = 0; i < n; i++) {
            tempEigenvalues[i] = koerper.getZero();
        }

        for (int j = 0; j < n; j++) {
            int its = 0;
            int m;
            do {
                for (m = j; m < n - 1; m++) {
                    Dfp delta = (tempEigenvalues[m].abs())
                            .add(tempEigenvalues[m + 1].abs());
                    if ((e[m].abs()).add(delta).equals(delta)) {
                        break;
                    }
                }
                if (m != j) {
                    if (its == maxIter) {
                        throw new MaxCountExceededException(
                                LocalizedFormats.CONVERGENCE_FAILED, maxIter);
                    }
                    its++;
                    Dfp q = (tempEigenvalues[j + 1]
                            .subtract(tempEigenvalues[j]))
                            .divide(e[j].multiply(2));
                    Dfp t = (koerper.getOne().add(q.pow(2))).sqrt();
                    if (q.lessThan(koerper.getZero())) {
                        q = tempEigenvalues[m].subtract(tempEigenvalues[j])
                                .add(e[j].divide(q.subtract(t)));
                    } else {
                        q = tempEigenvalues[m].subtract(tempEigenvalues[j])
                                .add(e[j].divide(q.add(t)));
                    }
                    Dfp u = koerper.getZero();
                    Dfp s = koerper.getOne();
                    Dfp c = koerper.getOne();
                    int i;
                    for (i = m - 1; i >= j; i--) {
                        Dfp p = s.multiply(e[i]);
                        Dfp h = c.multiply(e[i]);
                        if ((p).abs().greaterThan(q.abs()) || 
                                (p).abs().equals(q.abs())) {
                            c = q.divide(p);
                            t = (c.pow(2).add(koerper.getOne())).sqrt();
                            e[i + 1] = p.multiply(t);
                            s = t.reciprocal();
                            c = c.multiply(s);
                        } else {
                            s = p.divide(q);
                            t = (s.pow(2).add(koerper.getOne())).sqrt();
                            e[i + 1] = q.multiply(t);
                            c = t.reciprocal();
                            s = s.multiply(c);
                        }
                        if (e[i + 1].isZero()) {
                            tempEigenvalues[i + 1] = tempEigenvalues[i + 1]
                                    .subtract(u);
                            e[m] = koerper.getZero();
                            break;
                        }
                        q = tempEigenvalues[i + 1].subtract(u);
                        t = ((tempEigenvalues[i].subtract(q)).multiply(s))
                                .add(c.multiply(h).multiply(2));
                        u = s.multiply(t);
                        tempEigenvalues[i + 1] = q.add(u);
                        q = (c.multiply(t)).subtract(h);
                    }
                    if (t.isZero() && i >= j) {
                        continue;
                    }
                    tempEigenvalues[j] = tempEigenvalues[j].subtract(u);
                    e[j] = q;
                    e[m] = koerper.getZero();
                }
            } while (m != j);
        }
        /* Sort the eigen values in increase order. */
        for (int i = 0; i < n; i++) {
            int k = i;
            Dfp p = tempEigenvalues[i];
            for (int j = i + 1; j < n; j++) {
                if (tempEigenvalues[j].lessThan(p)) {
                    k = j;
                    p = tempEigenvalues[j];
                }
            }
            if (k != i) {
                tempEigenvalues[k] = tempEigenvalues[i];
                tempEigenvalues[i] = p;
            }
        }
        return tempEigenvalues;
    }

    /**
     * Gets a copy of the eigenvalues of the original matrix.
     * @return a copy of the eigenvalues of the original matrix.
     */
    public Dfp[] getEigenvalues() {
        return eigenvalues.clone();
    }
}
