package de.charite.compbio.exomiser.core.prioritisers;

/**
 * BOQA score in HPO
 *
 * @author Sebastian Koehler
 * @version 0.01 (20 April, 2016).
 */
public class Phenix2PriorityResult implements PriorityResult {

    /**
     * The semantic similarity score as implemented in PhenIX (also know as
     * Phenomizer). Note that this is not the p-value methodology in that paper,
     * but merely the simple semantic similarity score.
     */
    private double boqaProbability;

    private static double NORMALIZATION_FACTOR = 1f;

    public static void setNormalizationFactor(double factor) {
        NORMALIZATION_FACTOR = factor;
    }

    public Phenix2PriorityResult(double boqaProbability) {
        this.boqaProbability = boqaProbability;
    }

    @Override
    public PriorityType getPriorityType() {
        return PriorityType.PHENIX2_PRIORITY;
    }

    /**
     * @see exomizer.priority.IRelevanceScore#getRelevanceScore
     * @return the HPO semantic similarity score calculated via Phenomizer.
     */
    @Override
    public float getScore() {
        return (float) (boqaProbability * NORMALIZATION_FACTOR);
    }

    /**
     * @see exomizer.filter.Triage#getHTMLCode()
     */
    @Override
    public String getHTMLCode() {
        return String.format("<dl><dt>PhenIX2 BOQA score: %.2f</dt></dl>",
                this.boqaProbability);
    }

}
