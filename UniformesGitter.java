import org.apache.commons.math3.dfp.Dfp;

/**
 * Erzeugt ein Feld, das für einen übergebenen Paramter $l$
 * äquidistante Punkte $\xi_i, i = 0, \hdots, l$, im Intervall $[s, t]$
 * zurückgibt, mit $\xi_0 = s, \xi_l = t$.
 */
public class UniformesGitter implements Gitter<Dfp> {

    /**
     * Das Feld {@code Dfp[] xi} enthält die $l + 1$ Punkte.
     */
    final Dfp[] xi;

    /**
     * Erzeugt eine Instanz für ein $l > 0$, die alle Gitterknoten
     * $\xi_0, \hdots, \xi_l$ repräsentiert, also $l+1$ Knoten.
     *
     * @param l Index des letzten Gitterknotens.
     * @param s $= \xi_0$.
     * @param t $= \xi_l$.
     */
    public UniformesGitter(int l, Dfp s, Dfp t) {
        xi = setzePunkte(l, s, t);
    }

    /**
     * Erstellt ein {@code Dfp[]}-Array, das die
     * $\xi_i, i = 0, ..., l \subset [s, t]$ enthält.
     *
     * @param l Index des letzten Gitterknotens.
     * @param s $= \xi_0$.
     * @param t $= \xi_l$.
     * @return das {@code Dfp[]}-Array mit den $\xi_i$.
     */
    private Dfp[] setzePunkte(int l, Dfp s, Dfp t) {
        final Dfp c = t.subtract(s).divide(l);
        Dfp[] tempXi = new Dfp[l + 1];
        tempXi[0] = s;
        for (int j = 1; j <= l; j++) {
            tempXi[j] = tempXi[j - 1].add(c);
        }
        return tempXi;
    }

    /**
     * Gibt den $i$-ten Gitterknoten zurück, $i = 0, \hdots, l$.
     *
     * @return {@code xi[i]}
     */
    public Dfp getXi(int i) {
        return xi[0].newInstance(xi[i]);
    }

    /**
     * Gibt eine Kopie des Feldes {@code Dfp[] xi} der Gitterknoten
     * im Intervall $[s, t]$ zurück.
     *
     * @return eine Kopie von {@code Dfp[] xi}
     */
    public Dfp[] getXi() {
        return xi.clone();
    }

    /**
     * Gibt $l$, den Index des letzten Gitterknotens, zurück.
     *
     * @return $l$ Index des letzten Gitterknotens.
     */
    public int getL() {
        return xi.length;
    }
}
