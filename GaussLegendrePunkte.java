import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/**
 * Erzeugt ein Feld, das für einen übergebenen Parameter $k$ die
 * Nullstellen des $k$-ten Legendre-Polynoms bereitstellt.
 */
public class GaussLegendrePunkte {

    /**
     * Das Feld {@code Dfp[] rho} enthält die $k$ Punkte.
     */
    final Dfp[] rho;
    /**
     * Der Körper in dem die Gauß-Legendre-Punkte erzeugt werden sollen.
     */
    private final DfpField koerper;

    /**
     * Erzeugt eine Instanz für ein $k > 0$.
     *
     * @param k       die gewünschte Anzahl der Gauss-Legendre-Punkte.
     * @param koerper in dem die Punkte definiert sein sollen.
     */
    public GaussLegendrePunkte(int k, DfpField koerper) {
        this.koerper = koerper;
        rho = setzePunkte(k);
    }

    /**
     * Erstellt ein {@code double[]}-Array, das die
     * $\rho_j, j = 1, ..., k \subset (-1, 1)$ enthält.
     *
     * @param k Anzahl der Kollokationspunkte.
     * @return das {@code double[]}-Array mit den $\rho_j$.
     */
    private Dfp[] setzePunkte(int k) {
        Dfp[] neben = new Dfp[k];
        for (int j = 1; j <= k; j++) {
            neben[j - 1] = (koerper.getOne().negate().add(4 *
                    Math.pow(j, 2))).sqrt().reciprocal().multiply(j);
        }
        return new FieldEigenDecomposition(neben).getEigenvalues();
    }

    /**
     * Gibt den $j$-ten Kollokationspunkt zurück.
     *
     * @return $rho[j]$
     */
    public Dfp getRho(int j) {
        return rho[j];
    }

    /**
     * Gibt eine Kopie des Feldes {@code double[] rho} der Kollokationspunkte
     * im Intervall $(-1, 1)$ zurück.
     *
     * @return eine Kopie von {@code double[] rho}
     */
    public Dfp[] getRho() {
        return rho.clone();
    }

    /**
     * Gibt die Anzahl $k$ der erzeugten Kollokationspunkte zurück.
     *
     * @return $k$ Anzahl der Kollokationspunkte.
     */
    public int getK() {
        return rho.length;
    }
}
