import org.apache.commons.math3.analysis.RealFieldUnivariateFunction;
import org.apache.commons.math3.dfp.Dfp;

/**
 * Repräsentiert ein $\textsc{Bakhvalov}$-Gitter auf einem beliebigen Intervall $[s, t]$ mit
 * $l$ Teilintervallen in Form eines Feldes $\verb!Dfp []!$, welches die Knoten enthält.
 */
public class BakhvalovGitter implements Gitter<Dfp> {

    /**
     * Enthält die Knoten.
     */
    private final Dfp[] xiB;

    /**
     * Erzeugt ein $\textsc{Bakhvalov}$-Gitter $\xi^S = (\xi_i^S)$ mit $l+1$ Knoten für eine
     * Konvektionsdiffusionsgleichung.
     *
     * @param l       die Anzahl der Teilintervalle.
     * @param s       der erste Gitterknoten $\xi_0^S$.
     * @param t       der letzte Gitterknoten $\xi_l^S$.
     * @param q       ein Gitterparameter beschreibt größenordnungsmäßig den Anteil
     *                der Gitterknoten, welcher in der Grenzschicht liegt.
     * @param sigma   ein Gitterparameter beschreibt die Auflösung der
     *                Grenzschicht und wird typischerweise nahe der formalen
     *                Konvergenzordnung der verwendeten Methode gewählt.
     * @param beta    gegeben durch die Koeffizienten der Differentialgleichung
     *                $-\varepsilon u'' - bu' + cu = f$, für die das Gitter erzeugt wird, durch $b \geq \beta > 0$.
     * @param epsilon singulärer Störungsparameter für welchen das Gitter
     *                erzeugt werden soll.
     */
    public BakhvalovGitter(int l, final Dfp s, final Dfp t, final Dfp q,
                           Dfp sigma, Dfp beta, Dfp epsilon) {
        xiB = new Dfp[l + 1];
        xiB[0] = s;
        /* Wiederholt auftretender Ausdruck */
        final Dfp sigmaEpsilonDivBeta = sigma.multiply(epsilon).divide(beta);
        /*
         * Die gittererzeugende Funktion $\chi (r) := \frac{-\sigma \varepsilon}{\beta} \ln \frac{q - r}{q}$ für $r = [0, \tau]$
         * innerhalb der Grenzschicht.
         */
        RealFieldUnivariateFunction<Dfp> chi = new
                RealFieldUnivariateFunction<>() {
                    private final Dfp sigEpsDivBet = sigmaEpsilonDivBeta;

                    public Dfp value(Dfp x) {
                        return sigEpsDivBet.multiply(q.subtract(x).divide(q)
                                .log()).negate();
                    }
                };
        /*
         * Bestimmt $\tau$ iterativ als Lösung der Gleichung $\chi'(\tau_{i+1}) = \frac{1 - \chi(\tau_i)}{1 - \tau_i}$, die
         * genau dann existiert, wenn $\sigma \varepsilon < \beta q$. Andernfalls ist $tau = 0$ und
         * das Gitter global uniform.
         */
        Dfp tau = s.getZero();
        if (sigma.multiply(epsilon).lessThan(beta.multiply(q))) {
            Dfp tau_i;
            do {
                tau_i = tau;
                tau = q.subtract(sigmaEpsilonDivBeta.multiply(q.getOne()
                        .subtract(tau_i)).divide(q.getOne().subtract(chi
                        .value(tau_i))));
            } while (!(tau.subtract(tau_i).isZero()));
        }
        Dfp TMinusS = t.subtract(s), chiTau = chi.value(tau);
        for (int i = 1; i <= l; i++) {
            Dfp r_i = s.newInstance(i).divide(l);
            /* Bestimmung der $\xi_i$ mittels $\chi$, falls $r_i < \tau$, ...*/
            if (r_i.lessThan(tau)) {
                xiB[i] = s.add(chi.value(r_i).multiply(TMinusS));
            } else {
                /*
                 * ... und die Berechnung des ersten Gitterknotens außerhalb
                 * des Intervalls $[0, \tau]$ gefolgt von der additiven Bestimmung
                 * aller folgenden Knoten des uniformen Teilgitters.
                 */
                xiB[i] = s.add(chiTau.add(sigmaEpsilonDivBeta
                        .divide(q.subtract(tau)).multiply(r_i
                                .subtract(tau))).multiply(TMinusS));
                Dfp temp = TMinusS.subtract(xiB[i]).divide(l - i);
                for (int j = i + 1; j <= l; j++) {
                    xiB[j] = xiB[j - 1].add(temp);
                }
                break;
            }
        }
    }

    /**
     * Gibt den $i$-ten Gitterknoten zurück, $i = 0, \hdots, l$.
     *
     * @return $\verb!xi[i]!$
     */
    public Dfp getXi(int i) {
        return xiB[0].newInstance(xiB[i]);
    }

    /**
     * Gibt eine Kopie des Feldes $\verb!Dfp[] xiB!$ der Gitterknoten
     * im Intervall $[s, t]$ zurück.
     *
     * @return eine Kopie von $\verb!Dfp[] xiB!$
     */
    public Dfp[] getXi() {
        return xiB.clone();
    }
}
