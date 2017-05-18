import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;

/**
 * Erzeugt ein Feld, das für einen übergebenen Paramter $k$ 
 * äquidistante Punkte im Intervall $(-1, 1)$ zurückgibt.
 */
public class KollokationsPunkte {
    
    /**
     * Der Körper in dem die Punkte erzeugt werden soll.
     */
    private final DfpField koerper;
   /**
    * Das Feld {@code double[] rho} enthält die $k$ Punkte.
    */
    final Dfp[] rho;

    /**
     * Erzeugt eine Instanz für ein $k > 0$.
     * @param k
     */
    public KollokationsPunkte(int k, DfpField koerper) {
        this.koerper = koerper;
        rho = setzePunkte(k);
    }

  /**
   * Erstellt ein {@code double[]}-Array, dass die
   * $\rho_j, j = 1, ..., k \subset (-1, 1)$ enthält, aus denen dann
   * die $\tau_j$ berechnet werden können.
   * @param k Anzahl der Kollokationspunkte.
   * @return das {@code double[]}-Array mit den $\rho_j$.
   */
    private Dfp[] setzePunkte(int k) {
        Dfp[]  tempRho = new Dfp[k];
        Dfp temp = koerper.getOne().add(k).divide(koerper.getTwo());
        for (int j = 1; j <= tempRho.length; j++) {
            tempRho[j-1] = temp.reciprocal().multiply(j).subtract(1);
        };
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
     * Testprozedur für die Korrektheit der Implementierung
     * @param args
     */
    public static void main (String[] args) {        
        /**
         * Testen zweier Beispiele in verschiedenen Modi durch Ausgabe
         * interessierender Daten auf der Konsole.
         */
        final int n = 12;
        KollokationsPunkte rho;
        for (int k = 1; k <= n; k++) {
            System.out.println("k = " + k + "\n");
            rho = new KollokationsPunkte(k, new DfpField(30));
            for (int j = 1; j <= k; j++) {
                System.out.println(rho.getRho(j-1));
            }
            System.out.println();
        }
    }
}
