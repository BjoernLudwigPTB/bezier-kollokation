import org.apache.commons.math3.dfp.Dfp;

/**
 * Repräsentiert ein $\textsc{Shishkin}$-Gitter auf einem beliebigen Intervall $[s, t]$ mit
 * $l$ Teilintervallen in Form eines $\verb!Dfp []!$ Feldes, welches die Knoten enthält.
 */
public class ShishkinGitter implements Gitter<Dfp> {

    /**
     * Enthält die Knoten.
     */
    private final Dfp[] xiS;

    /**
     * Erzeugt ein $\textsc{Shishkin}$-Gitter $\xi^S = (\xi_i^S)$ mit $l+1$ Knoten für eine
     * Konvektionsdiffusionsgleichung.
     *
     * @param l       die Anzahl der Gitterintervalle.
     * @param s       der erste Gitterknoten $\xi_0^S$.
     * @param t       der letzte Gitterknoten $\xi_l^S$.
     * @param q       ein Gitterparameter beschreibt größenordnungsmäßig den
     *                Anteil der Gitterknoten, welcher in der Grenzschicht liegt.
     * @param sigma   ein Gitterparameter beschreibt die Auflösung der
     *                Grenzschicht und wird typischerweise nahe der formalen
     *                Konvergenzordnung der verwendeten Methode gewählt.
     * @param beta    gegeben durch die Koeffizienten der Differentialgleichung
     *                $-\varepsilon u'' - bu' + cu = f$, für die das Gitter erzeugt wird, durch $b \geq \beta > 0$.
     * @param epsilon singulärer Störungsparameter für welchen das Gitter
     *                erzeugt werden soll.
     */
    public ShishkinGitter(int l, Dfp s, Dfp t, Dfp q, Dfp sigma, Dfp beta,
                          Dfp epsilon) {
        xiS = new Dfp[l + 1];
        xiS[0] = s;
        int qL = q.newInstance(q.multiply(l).floor()).intValue();
        /* Berechnung des Gitterübergangpunktes $\tau$. */
        Dfp tau = sigma.multiply(epsilon).divide(beta).multiply(Math.log(l));
        if (q.lessThan(tau)) tau = q;
        /* Berechnung der lokal uniformen Teilgitter... */
        Dfp TMinusS = t.subtract(s);
        {
            /*... im Intervall $[s, s + \tau (t - s)]$. */
            Dfp tauTMinusSDivQL = tau.multiply(TMinusS).divide(qL);
            for (int i = 1; i <= qL; i++) {
                xiS[i] = xiS[i - 1].add(tauTMinusSDivQL);
            }
        }
        {
            /*... im Intervall $[s + \tau (t - s), t], t]$. */
            Dfp einsMinustauTMinusSDivQL = tau.getOne().subtract(tau)
                    .multiply(TMinusS).divide(l - qL);
            for (int i = qL + 1; i <= l; i++) {
                xiS[i] = xiS[i - 1].add(einsMinustauTMinusSDivQL);
            }
        }
    }

    /**
     * Erzeugt ein $\textsc{Shishkin}$-Gitter $\xi^S = (\xi_i^S)$ mit $l+1$ Knoten für eine
     * Reaktionsdiffusionsgleichung.
     *
     * @param l       die Anzahl der Teilintervalle.
     * @param s       der erste Gitterknoten $\xi_0^S$.
     * @param t       der letzte Gitterknoten $\xi_l^S$.
     * @param q_0     ein Gitterparameter beschreibt größenordnungsmäßig den
     *                Anteil der Gitterknoten, welcher in der Grenzschicht bei $x=s$ liegt.
     * @param q_1     ein Gitterparameter beschreibt größenordnungsmäßig den
     *                Anteil der Gitterknoten, welcher in der Grenzschicht bei $x=t$ liegt.
     * @param sigma_0 ein Gitterparameter beschreibt die Auflösung der
     *                Grenzschicht bei $x=s$ und wird typischerweise nahe der formalen
     *                Konvergenzordnung der verwendeten Methode gewählt.
     * @param sigma_1 ein Gitterparameter beschreibt die Auflösung der
     *                Grenzschicht bei $x=t$ und wird typischerweise nahe der formalen
     *                Konvergenzordnung der verwendeten Methode gewählt.
     * @param gamma   gegeben durch die Koeffizienten der Differentialgleichung
     *                $-\varepsilon u'' + cu = f$, für die das Gitter erzeugt wird, durch $c \geq \gamma^2, \gamma > 0$.
     * @param epsilon singulärer Störungsparameter für welchen das Gitter
     *                erzeugt werden soll.
     */
    public ShishkinGitter(int l, Dfp s, Dfp t, Dfp q_0, Dfp q_1, Dfp sigma_0,
                          Dfp sigma_1, Dfp gamma, Dfp epsilon) {
        xiS = new Dfp[l + 1];
        xiS[0] = s;
        int q_0L = q_0.newInstance(q_0.multiply(l).floor()).intValue(),
                q_1L = q_1.newInstance(q_1.multiply(l).floor()).intValue();
        Dfp tau_0, tau_1;
        /* Berechnung der Gitterübergangpunkte $\tau_0$ und $\tau_1$. */
        {
            Dfp temp = epsilon.divide(gamma).multiply(Math.log(l));
            tau_0 = sigma_0.multiply(temp);
            tau_1 = sigma_1.multiply(temp);
        }
        if (q_0.lessThan(tau_0)) tau_0 = q_0;
        if (q_1.lessThan(tau_1)) tau_1 = q_1;
        /* Berechnung der lokal uniformen Teilgitter... */
        Dfp TMinusS = t.subtract(s);
        {
            /*... im Intervall $[s, s + \tau_0 (t - s)]$. */
            Dfp tau_0TMinusSDivQ_0L = tau_0.multiply(TMinusS).divide(q_0L);
            for (int i = 1; i <= q_0L; i++) {
                xiS[i] = xiS[i - 1].add(tau_0TMinusSDivQ_0L);
            }
        }
        {
            /*... im Intervall $[s + \tau_0 (t - s), t - \tau_1 (t - s)]$. */
            Dfp einsMinustau_0Plustau1TMinusSDivQL = tau_0.getOne().subtract(
                    tau_0).subtract(tau_1).multiply(TMinusS)
                    .divide(l - (q_0L + q_1L));
            for (int i = q_0L + 1; i <= l - q_1L; i++) {
                xiS[i] = xiS[i - 1].add(einsMinustau_0Plustau1TMinusSDivQL);
            }
        }
        {
            /*... im Intervall $[t - \tau_1 (t - s)], t]$. */
            System.out.println();
            Dfp tau_1TMinusSDivQ_1L = tau_1.multiply(TMinusS).divide(q_1L);
            for (int i = l - q_1L + 1; i <= l; i++) {
                xiS[i] = xiS[i - 1].add(tau_1TMinusSDivQ_1L);
            }
        }
    }

    /**
     * Gibt den $i$-ten Gitterknoten zurück, $i = 0, \hdots, l$.
     *
     * @return $\verb!xi[i]!$
     */
    public Dfp getXi(int i) {
        return xiS[0].newInstance(xiS[i]);
    }

    /**
     * Gibt eine Kopie des Feldes $\verb!Dfp[] xiS!$ der Gitterknoten zurück.
     *
     * @return eine Kopie von $\verb!Dfp[] xiS!$
     */
    public Dfp[] getXi() {
        return xiS.clone();
    }
}
