import org.apache.commons.math3.dfp.Dfp;
import org.apache.commons.math3.dfp.DfpField;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.*;

/**
 * Implementierung des Algorithmus $\verb!bandet1!$ von Martin und Wilkinson (
 * Numer. Math. 9 (1967) pp. 279--301) zur direkten Lösung eines linearen
 * Gleichungssystems mit einer asymmetrischen Bandmatrix in der Version von
 * de Boor (A practical guide to splines (1978) pp. 370--372,
 * $\verb!SUBROUTINE CWIDTH!$), der den Algorithmus noch auf die spezielle
 * fast blockdiagonale Struktur der bei Spline-Kollokationen auftretenden
 * Matrizen angepasst hat. Lediglich die skalierte Spalten-Pivotsuche nehmen
 * wir bezogen auf die Zeilensummennorm vor und nicht, wie im Original, in
 * Bezug auf die Maximumnorm, weil diese als eine von zwei Pivotstrategien,
 * welche die Restmatrix einbeziehen, in $\textsc{T. Linß}$ (Numerische Mathematik I
 * (2015) pp. 77-78) genannt wird. Die Blockstruktur der Kollokationsmatrix
 * ist dabei auf eine Kollokationsmatrix zur Bestimmung eines $C^1$-Splines
 * ausgelegt, der eine Zwei-Punkt-Randwertaufgabe löst. Das heißt die Matrix 
 * $A$ beginnt mit einer Zeile zur Randbedingungen am Anfang des betrachteten
 * Intervalls, gefolgt von $k$ Kollokationsbedingungen, auf die wiederum bis
 * zum vorletzten der $l-1$ Blöcke $2$ Stetigkeitsbedingungen folgen und an
 * die sich im letzten Block die Randbedingung am rechten Ende des
 * betrachteten Intervalls anschließt.
 */
public class FieldBlockDecomposition {

    /** Körper auf dem die Zerlegung definiert ist. */
    private final DfpField koerper;
    /** Matrix $A$ im Gleichungssystem $Ax = b$. */
    private final Dfp[][] a;
    /** Anzahl der Blöcke in der Matrix $A$. */
    private final int k;
    /** 
     * Die Anzahl der Kollokationsbedingungen je Block. */
    private final int l;
    
    /**
     * Erzeugt eine Instanz des Zerlegers der Matrix $A$.
     * @param a die Matrix $A$ für die $Ax = b$ gelöst werden soll.
     * @param b der Vektor $b$ für den $Ax = b$ gelöst werden soll.
     * @param k die Anzahl der Kollokationsbedingungen mit der die Matrix
     * erstellt wurde.
     * @param l die Anzahl der Gitterintervalle. 
     */
    public FieldBlockDecomposition (Dfp[][] a, final int k,
            final int l) {
        if (a.length != l * (k + 2) && a[0].length == k+2)
            throw new DimensionMismatchException(a.length, l * (k+2));
        if (k == 1 && l > 1) {
            if (a[0].length != 4)
                throw new DimensionMismatchException(a[0].length, k+2);
        } else {
            if (a[0].length != k+2)
                throw new DimensionMismatchException(a[0].length, k+2);
        }
        this.k = k;
        this.l = l;
        koerper = (DfpField) a[0][0].getField();
        this.a = a;
    }
    
    /**
     * Löst für einen übergebenen Vektor $b$ das Gleichungssystem $Ax = b$.
     * @param b die rechte Seite des Gleichungssystems.
     * @return Die Lösung $x$ des Gleichungssystems $Ax = b$.
     */
    public Dfp[] solve(Dfp[] b) {
        /*
         * Folgende Variablenbezeichnungen verwenden wir im
         * Vergleich zu $\textsc{de Boor}$ (1978):
         * W $=: a$;
         * B $=: b$;
         * NEQU $=: a.$length;
         * COLMAX $=:$ maxDominanz;
         * TEMP $=:$ maxDominanz;
         * AWILOD $=:$ tempDominanz;
         * ISTAR $=:$ indexMaxDominanz;
         * NCOLS $=: a[0].$length;
         * LASTCL $=:$ letzteSpalteBlock;
         * NROWAD $=:$ structure[i][0];
         * NEXTEQ $=:$ letzteGleichungBlock + j;
         * INTEGS $=:$ structure;
         * NBLOKS $=: 2 * l - 1 = \operatorname{structure.length} - 1$;
         * D $=:$ zeilenSumme;
         * X $=:$ tempX;
         * IFLAG nicht vorhanden, weil wir keine Determinante berechnen;
         * IPVTEQ $=:$ pivot;
         * LASTEQ $=:$ letzteGleichungBlock;
         * PVTL $=:$ pivot $+1$;
         * I $=: i$;
         * II $=:$ restBlock;
         * jMax $=:$ restBlock;
         * ICOUNT $=: j$;
         * RATIO $=:$ ratio;
         * SUM $=:$ tempSumme;
         */
        if (b.length != a.length)
            throw new DimensionMismatchException(b.length, a.length);
        Dfp[] tempX = new Dfp[b.length];
        Dfp[] zeilenSumme = new Dfp[b.length];
        /*
         * Nimmt die Struktur der Blockmatrix auf. Dabei ist $\verb!structure[i][0]!$
         * die Anzahl der Zeilen in Block $i$ und $\verb!structure[i][1]!$ die Anzahl
         * der möglichen Pivotschritte in Block $i$, bevor möglichweise eine
         * Zeile aus dem nächsten Block in den Block hineinpivotiert wird.
         */
        int[][] structure;
        if (l == 1) {
            /*
             * Im Falle, dass nur ein Block vorhanden ist, wird über alle
             * Zeilen pivotiert.
             */
            structure = new int[][] {{k+2, k+2}};
        } else {
            /*
             * Sonst bilden die jeweils $k$ Kollokationsbedingungen jedes
             * Gitterintervalls einen Block, am oberen und unteren Rand der
             * Matrix zusammen mit den Randbedingungen, und jeweils die
             * zwei Zeilen der Stetigkeitsbedingungen zwischen den
             * Kollokationsbedingungen.
             */
            structure = new int[2 * l - 1][2];
            structure[0] = new int[] {k+1, k};
            for (int i = 1; i < l; i++) {
                structure[2*i-1] = new int[] {2, 2};
                structure[2*i] = new int[] {k, k};
            }
            /*
             * Für $k = 1$ entsteht bei der blockweisen Bearbeitung die
             * Besonderheit, dass die Blöcke der Länge $k + 2 = 3$ nicht mehr 
             * die Stetigkeitsbedingungen aufnehmen können, da sich diese
             * jeweils auf vier Bézierpunkte beziehen. Deshalb wird in diesem
             * Fall die Matrix geringfügig umsortiert und $\verb!structure!$
             * angepasst. Der vorletzte Block besteht dann lediglich aus der
             * $C^1$-Stetigkeitsbedingung und die $C^0$_Stetigkeitsbedingung
             * rutscht in den letzten Block. Dafür muss dann die zugehörige
             * Zeile in der Matrix um eins nach links und die letzte
             * Kollokationsbedingung eins nach rechts verschoben werden.
             */
            if (k == 1) {
                structure[structure.length-2] = structure[structure.length-1];
                structure[structure.length-1] = new int[] {3, 4};
                a[(l - 1) * 3][0] = a[(l - 1) * 3][1];
                a[(l - 1) * 3][1] = a[(l - 1) * 3][2];
                a[(l - 1) * 3][2] = koerper.getZero();
                a[(l - 1) * 3 + 1][3] = a[(l - 1) * 3 + 1][2];
                a[(l - 1) * 3 + 1][2] = a[(l - 1) * 3 + 1][1];
                a[(l - 1) * 3 + 1][1] = a[(l - 1) * 3 + 1][0];
                a[(l - 1) * 3 + 1][0] = koerper.getZero();
            } else
                structure[structure.length-1] = new int[] {k+1, k+2};
        }       
        int pivot = -1, letzteGleichungBlock = 0;
        {
            for (int i = 0; i < structure.length; i++) {
                { 
                    /*
                     * In diesem Bereich wird für die skalierte Spalten-
                     * Pivot-Suche für jeden in die Elimination erstmals 
                     * eintretenden Block jeweils für alle Zeilen die Summe
                     * der Beträge ihrer Einträge berechnet.
                     */
                    for (int j = 0; j < structure[i][0]; j++) {
                        Dfp tempZeilenSumme = koerper.getZero();
                        for (int m = 0; m < a[i].length; m++) {
                            tempZeilenSumme = tempZeilenSumme
                                    .add(a[letzteGleichungBlock+ j][m].abs());
                        }
                        if (tempZeilenSumme == koerper.getZero())
                            throw new SingularMatrixException();
                        zeilenSumme[letzteGleichungBlock+j] = tempZeilenSumme;
                    } 
                }
                /* Bezieht den nächsten Block in die Betrachtung ein. */
                letzteGleichungBlock += structure[i][0];
                /* 
                 * Im aktuellen Block wird jeweils mit voller Spaltenanzahl
                 * gestartet und pro Pivotschritt eine Spalte weniger
                 * berücksichtigt. Die potentiellen Nichtnullelemente des
                 * Blocks stehen dann immer beginnend mit Spalte eins in $a$.
                 */
                int letzteSpalteBlock = a[0].length;
                for (int j = 1; j <= structure[i][1]; j++) {
                    /*
                     * Wandert durch den aktuellen Block bis möglicherweise
                     * Zeilen aus dem nächsten Block hineinpivotiert werden.
                     */
                    pivot++;
                    if (pivot < letzteGleichungBlock) {
                        {
                            /*
                             * In diesem Bereich wird die skalierte Spalten-
                             * Pivot-Suche durchgeführt und gegebenenfalls der
                             * Austausch der Zeilen vorgenommen.
                             */
                            int indexMaxDominanz = pivot;
                            Dfp maxDominanz = a[pivot][0].abs()
                                    .divide(zeilenSumme[pivot]);
                            for (int m = pivot + 1; m < letzteGleichungBlock;
                                    m++) {
                                Dfp tempDominanz = a[m][0].abs()
                                        .divide(zeilenSumme[m]);
                                if (tempDominanz.greaterThan(maxDominanz)) {
                                    maxDominanz = tempDominanz;
                                    indexMaxDominanz = m;
                                }
                            }
                            if (indexMaxDominanz != pivot) {
                                maxDominanz = zeilenSumme[indexMaxDominanz];
                                zeilenSumme[indexMaxDominanz] = 
                                        zeilenSumme[pivot];
                                zeilenSumme[pivot] = maxDominanz;
                                maxDominanz = b[indexMaxDominanz];
                                b[indexMaxDominanz] = b[pivot];
                                b[pivot] = maxDominanz;
                                for (int m = 0; m < letzteSpalteBlock; m++) {
                                    maxDominanz = a[indexMaxDominanz][m];
                                    a[indexMaxDominanz][m] = a[pivot][m];
                                    a[pivot][m] = maxDominanz;
                                }
                            }
                        }
                        /*
                         * Durchführung der Elimination und verschieben aller
                         * verbleibenden Elemente um eine Spalte nach links.
                         */
                        for (int restBlock = pivot + 1;
                                restBlock<letzteGleichungBlock; restBlock++) {
                            Dfp ratio = a[restBlock][0].divide(a[pivot][0]);
                            for (int m = 1; m < letzteSpalteBlock; m++) {
                                a[restBlock][m-1] = a[restBlock][m]
                                        .subtract(ratio
                                                .multiply(a[pivot][m]));
                            }
                            a[restBlock][letzteSpalteBlock-1] = koerper
                                    .getZero();
                            b[restBlock] = b[restBlock]
                                    .subtract(ratio.multiply(b[pivot]));
                        }
                        letzteSpalteBlock--;
                    } else if (a[pivot][0].isZero())
                        throw new SingularMatrixException();
                }    
            }
        }
        /* Initialisierung des Lösungsvektors. */
        for (int i = 0; i < tempX.length; i++) {
            tempX[i] = koerper.getZero();
        }
        /* Rückwärtsersetzung unter Berücksichtigung der Struktur von $A$. */
        for (int i = structure.length - 1; i >= 0; i--) {
            int restBlock = a[0].length - structure[i][1];
            for (int j = 0; j < structure[i][1]; j++) {
                Dfp tempSumme = koerper.getZero();
                for (int m = 1; m <= restBlock; m++) {
                    tempSumme = tempSumme.add(tempX[pivot+m]
                            .multiply(a[pivot][m]));
                }
                tempX[pivot]=b[pivot].subtract(tempSumme)
                        .divide(a[pivot][0]);
                restBlock++;
                pivot--;
            }
        }
        return tempX;
    }
}
