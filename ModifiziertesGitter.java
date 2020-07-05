import org.apache.commons.math3.dfp.Dfp;

public class ModifiziertesGitter implements Gitter<Dfp> {

    /**
     * Enthält die Gitterknoten.
     */
    private final Dfp[] xi;

    /**
     * Erzeugt zu einem übergebenen Gitter ein Gitter, welches $r$-mal feiner
     * aufgelöst ist, in dem jedes Gitterintervall in $r$  äquidistante
     * Gitterintervalle zerlegt wird.
     *
     * @param gitter das feiner aufgelöst werden soll.
     * @param r      mal feiner wird das neue Gitter aufgelöst.
     */
    public ModifiziertesGitter(Gitter<Dfp> gitter, int r) {
        Dfp[] tempXi = gitter.getXi();
        int tempL = r * (tempXi.length - 1) + 1;
        xi = new Dfp[tempL];
        xi[0] = tempXi[0];
        for (int i = 1; i < tempXi.length; i++) {
            xi[r * i] = gitter.getXi(i);
            Dfp temp = tempXi[i].subtract(tempXi[i - 1]).divide(r);
            for (int j = 1; j < r; j++) {
                xi[r * (i - 1) + j] = xi[r * (i - 1) + (j - 1)].add(temp);
            }
        }
    }

    /* (non-Javadoc)
     * @see BakhvalovGitter#getXi(int)
     */
    @Override
    public Dfp getXi(int i) {
        return xi[0].newInstance(xi[i]);
    }

    /* (non-Javadoc)
     * @see BakhvalovGitter#getXi()
     */
    @Override
    public Dfp[] getXi() {
        return xi.clone();
    }
}
