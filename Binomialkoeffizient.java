/**
 * Eine Instanz enthält für einen Paramter $n$ alle Werte von
 * $n$ über $k, k = 0, ..., n$.
 */
public class Binomialkoeffizient {

    /** Enthält den Wert $n$. */
    private int n;
    /** Enthält die berechneten Werte. */
    private int[] binom;
    
    /**
     * Erzeugt eine Instanz, die alle entsprechenden Werte
     * $n$ über $k, k = 0, ..., n$ vorhält.
     */
    public Binomialkoeffizient(int n) {
        this.n = n;
        binom = new int[n+1];
        /**
         * Bei der Berechnung der Werte wird die Symmetrie des
         * Binomialkoeffizienten ausgenutzt und für $k > n - k$ der
         * bereits berechnete Wert für $n - k$ genutzt.
         */
        for (int k = 0; k < binom.length; k++) {
            if (k > n - k)
                binom[k] = binom[n - k];
            else
                binom[k] = berechneNUeberK(k);
        }
    }

    /**
     * Berechnet die Werte des Binomialkoeffizienten mit der Variablen
     * $n$ und dem übergebenen Paramter.
     * @param k für das $n$ über $k$ berechnet werden  soll.
     * @return $n$ über $k$
     */
    private int berechneNUeberK(int k) {
        int b = 1;
        for (int j = 1, m = n; j <= k; j++, m--)
            b = b * m/j;
        return b;
    }
    
    /**
     * Gibt den Wert $n$ über $k$ zurück.
     * @param k für das $n$ über $k$ zurückgegeben werden soll.
     * @return $n$ über $k$.
     */
    public int getUeber(int k) {
        return binom[k];
    }
    
    /**
     * Gibt alle Werte $n$ über $k, k = 0, ..., n$ zurück.
     * @param k für das $n$ über $k$ zurückgegeben werden soll.
     * @return {@code double[] (n über 0, ..., n über n})}.
     */
    public int[] getBinom() {
        return (int[]) binom.clone();
    }
    
    /**
     * Gibt den Wert zurück, für den die Werte des Binomialkoeffizienten
     * abgelegt sind.
     * @return $n$, für das alle $n$ über $k, k = 0, ..., n$ abgelegt
     * sind.
     */
    public int getN() {
        return n;
    }
}
