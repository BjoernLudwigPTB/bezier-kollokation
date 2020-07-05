import org.apache.commons.math3.RealFieldElement;

/**
 * Parametrisiertes Interface für die Definition von Gittern.
 */
public interface Gitter<T extends RealFieldElement<T>> {

    /**
     * Gibt den $i$-ten Gitterknoten zurück, $i = 0, \hdots, l$.
     *
     * @return $\verb!xi[i]!$
     */
    T getXi(int i);

    /**
     * Gibt eine Kopie des Feldes $\verb!Dfp[] xiS!$ der Gitterknoten
     * im Intervall $[s, t]$ zurück.
     *
     * @return eine Kopie von $\verb!Dfp[] xiS!$
     */
    T[] getXi();
}
