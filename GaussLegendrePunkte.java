import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.util.FastMath;

/**
 * Erzeugt ein Feld, das für einen übergebenen Paramter $k$ die
 * Nullstellen des $k$-ten Legendre-Polynoms bereitstellt.
 */
public class GaussLegendrePunkte {
    
    /** Der Körper in dem die Gauß-Legendre-Punkte erzeugt werden sollen. */
    private final DfpField koerper;
   /** Das Feld {@code Dfp[] rho} enthält die $k$ Punkte. */
    final Dfp[] rho;

    /**
     * Erzeugt eine Instanz für ein $k > 0$.
     * @param k
     */
    public GaussLegendrePunkte(int k, DfpField koerper) {
        this.koerper = koerper;
        rho = setzePunkte(k);
    }

  /**
   * Erstellt ein {@code double[]}-Array, das die
   * $\rho_j, j = 1, ..., k \subset (-1, 1)$ enthält, aus denen dann
   * die $\tau_j$ berechnet werden können.
   * @param k Anzahl der Kollokationspunkte.
   * @return das {@code double[]}-Array mit den $\rho_j$.
   */
    private Dfp[] setzePunkte(int k) {
        Dfp[]  tempRho, neben = new Dfp[k];
        for (int j = 1; j <= k; j++) {
            neben[j-1] = (koerper.getOne().negate().add(4 * 
                    FastMath.pow(j, 2))).sqrt().reciprocal().multiply(j);
        };
        tempRho = new FieldEigenDecomposition(neben).getEigenvalues();
        ArrayUtils.reverse(tempRho);
        return tempRho;
    }

    /**
     * Gibt den $j$-ten Kollokationspunkt zurück.
     * @return $rho[j]$
     */
    public Dfp getRho(int j) {
        return rho[j];
    }

    /**
     * Gibt eine Kopie des Feldes {@code double[] rho} der Kollokationspunkte
     * im Intervall $(-1, 1)$ zurück.
     * @return eine Kopie von {@code double[] rho}
     */
    public Dfp[] getRho() {
        return rho.clone();
    }

    /**
     * Gibt die Anzahl $k$ der erzeugten Kollokationspunkte zurück.
     * @return $k$ Anzahl der Kollokationspunkte.
     */
    public int getK() {
        return rho.length;
    }

    /**
     * Testprozedur für die Korrektheit der Implementierung.
     * @param args
     */
    public static void main (String[] args) {
        final DfpField koerper = new DfpField(50);
        final int n = 12;
        GaussLegendrePunkte rho1;
        for (int k = 1; k <= n; k++) {
            System.out.println("k = " + k + "\n");
            rho1 = new GaussLegendrePunkte(k, koerper);
            for (int j = 1; j <= k; j++) {
                System.out.println(rho1.getRho(j-1));
            }
            System.out.println();
        }
    }
}
