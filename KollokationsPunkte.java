import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/**
 * Erzeugt ein Feld, das für einen übergebenen Parameter $k$
 * äquidistante Punkte im Intervall $(-1, 1)$ zurück gibt.
 */
public class KollokationsPunkte {
    
    /** Der Körper in dem die Punkte erzeugt werden soll. */
    private final DfpField koerper;
   /** Das Feld {@code Dfp[] rho} enthält die $k$ Punkte. */
    final Dfp[] rho;

    /**
     * Erzeugt eine Instanz für ein $k > 0$.
     * @param k Anzahl der Kollokationspunkte.
     */
    public KollokationsPunkte(int k, DfpField koerper) {
        this.koerper = koerper;
        rho = setzePunkte(k);
    }

  /**
   * Erstellt ein {@code Dfp[]}-Array, das die
   * $\rho_j, j = 1, ..., k \subset (-1, 1)$ enthält.
   * @param k Anzahl der Kollokationspunkte.
   * @return das {@code Dfp[]}-Array mit den $\rho_j$.
   */
    private Dfp[] setzePunkte(int k) {
        Dfp[]  tempRho = new Dfp[k];
        Dfp temp = koerper.getOne().add(k).divide(koerper.getTwo());
        for (int j = 1; j <= k; j++) {
            tempRho[j-1] = temp.reciprocal().multiply(j).subtract(1);
        }
        return tempRho;
    }

    /**
     * Gibt den $j$-ten Kollokationspunkt zurück.
     * @return {@code rho[j]}
     */
    public Dfp getRho(int j) {
        return rho[j];
    }

    /**
     * Gibt eine Kopie des Feldes {@code Dfp[] rho} der Kollokationspunkte
     * im Intervall $(-1, 1)$ zurück.
     * @return eine Kopie von {@code Dfp[] rho}
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
}
