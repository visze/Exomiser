/*
 * The Exomiser - A tool to annotate and prioritize variants
 *
 * Copyright (C) 2012 - 2015  Charite Universitätsmedizin Berlin and Genome Research Ltd.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core.analysis.util;

import de.charite.compbio.exomiser.core.model.Gene;
import de.charite.compbio.exomiser.core.model.TopologicalDomain;
import de.charite.compbio.exomiser.core.model.VariantEvaluation;
import de.charite.compbio.exomiser.core.prioritisers.PriorityType;
import de.charite.compbio.jannovar.annotation.Annotation;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;

import static java.util.stream.Collectors.toList;

/**
 * Reassigns regulatory non-coding variants to the gene with the best phenotype score in a topological domain
 * (doi:10.1038/nature11082). 'Recent research shows that high-order chromosome structures make an important
 * contribution to enhancer functionality by triggering their physical interactions with target genes.'
 * (doi:10.1038/nature12753).
 *
 * @author Damian Smedley <damian.smedley@sanger.ac.uk>
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class GeneReassigner {

	private static final Logger logger = LoggerFactory.getLogger(GeneReassigner.class);

	private final PriorityType priorityType;
	private final ChromosomalRegionIndex<TopologicalDomain> tadIndex;
	private final Map<String, Gene> allGenes;

	/**
	 * @param priorityType
	 * @param tadIndex
	 */
	public GeneReassigner(PriorityType priorityType, Map<String, Gene> allGenes,
			ChromosomalRegionIndex<TopologicalDomain> tadIndex) {
		this.tadIndex = tadIndex;
		this.allGenes = allGenes;
		this.priorityType = priorityType;
		logger.info("Made new GeneReassigner for {}", priorityType);
	}

	public void reassignVariantToMostPhenotypicallySimilarGeneInTad(VariantEvaluation variantEvaluation) {
        switch (variantEvaluation.getVariantEffect()) {
        case CODING_TRANSCRIPT_INTRON_VARIANT:
        case CONSERVED_INTERGENIC_VARIANT:
        case CONSERVED_INTRON_VARIANT:
        case DOWNSTREAM_GENE_VARIANT:
        case INTERGENIC_REGION:
        case INTERGENIC_VARIANT:
        case INTRAGENIC_VARIANT:
        case INTRON_VARIANT:
        case NON_CODING_TRANSCRIPT_INTRON_VARIANT:
        case REGULATORY_REGION_VARIANT:
        case TF_BINDING_SITE_VARIANT:
        case UPSTREAM_GENE_VARIANT:
            logger.debug("Checking gene assignment for {} chr={} pos={}", variantEvaluation.getVariantEffect(), variantEvaluation.getChromosome(), variantEvaluation.getPosition());
            assignVariantToGeneWithHighestPhenotypeScore(variantEvaluation);
            break;
            default:
            	break;
        }
    }

	private void assignVariantToGeneWithHighestPhenotypeScore(VariantEvaluation variantEvaluation) {
		Gene geneWithHighestPhenotypeScore = null;
		Gene currentlyAssignedGene = allGenes.get(variantEvaluation.getGeneSymbol());
		// assign this to the variant's current gene as we don't necessarily want ALL the regulatory region variants to
		// clump into one gene.
		float bestScore = prioritiserScore(currentlyAssignedGene);
		List<String> genesInTad = getGenesInTadForVariant(variantEvaluation);

		for (String geneSymbol : genesInTad) {
			Gene gene = allGenes.get(geneSymbol);
			logger.debug("Checking gene {}", gene);
			if (gene != null && (gene.getPriorityResult(priorityType)) != null) {
				float geneScore = prioritiserScore(gene);
				logger.debug("Gene {} in TAD has score {}", geneSymbol, geneScore);
				if (geneScore > bestScore) {
					bestScore = geneScore;
					geneWithHighestPhenotypeScore = gene;
				}
			}
		}
		if (prioritiserScore(currentlyAssignedGene) == bestScore) {
			// don't move the assignment if there is nowhere better to go...
			return;
		}
		if (geneWithHighestPhenotypeScore == null) {
			return;
		}
		assignVariantToGene(variantEvaluation, geneWithHighestPhenotypeScore);
	}

	private float prioritiserScore(Gene gene) {
		// TODO: local variable, only one call to get? Otherwise will throw NPE if not guarded against.
		return (gene.getPriorityResult(priorityType) == null) ? 0f : gene.getPriorityResult(priorityType).getScore();
	}

	private List<String> getGenesInTadForVariant(VariantEvaluation variantEvaluation) {
		List<TopologicalDomain> tadsContainingVariant = tadIndex.getRegionsContainingVariant(variantEvaluation);

		return tadsContainingVariant.stream().map(TopologicalDomain::getGenes)
				.flatMap(geneMap -> geneMap.keySet().stream()).collect(toList());
	}

	private void assignVariantToGene(VariantEvaluation variantEvaluation, Gene gene) {
		logger.debug("Reassigning {} {}:{} {}->{} from {} to {}", variantEvaluation.getVariantEffect(),
				variantEvaluation.getChromosome(), variantEvaluation.getPosition(), variantEvaluation.getRef(),
				variantEvaluation.getAlt(), variantEvaluation.getGeneSymbol(), gene.getGeneSymbol());
		variantEvaluation.setEntrezGeneId(gene.getEntrezGeneID());
		variantEvaluation.setGeneSymbol(gene.getGeneSymbol());
		// given the physical ranges of topologically associated domains, the annotations are likely to be meaningless
		// once reassigned
		// TODO: check if the reassigned gene is actually in the annotation list and include this.
		variantEvaluation.setAnnotations(Collections.emptyList());
	}

	public void reassignGeneToMostPhenotypicallySimilarGeneInAnnotations(VariantEvaluation variantEvaluation) {
		List<Annotation> annotations = variantEvaluation.getAnnotations();
		Gene geneWithHighestPhenotypeScore = null;
		VariantEffect variantEffectForTopHit = null;
		Gene currentlyAssignedGene = allGenes.get(variantEvaluation.getGeneSymbol());
		float bestScore = 0f;
		if (currentlyAssignedGene == null) {
			// very rarely a variant just has a single annotation with no gene i.e. geneSymbol is .
			return;
		} else {
			bestScore = prioritiserScore(currentlyAssignedGene);
		}
		List<String> geneSymbols = new ArrayList<>();
		List<VariantEffect> variantEffects = new ArrayList<>();
		List<Annotation> newAnnotations = new ArrayList<>();
		for (Annotation a : annotations) {
			String geneSymbol = a.getGeneSymbol();
			geneSymbols.add(geneSymbol);
			variantEffects.add(a.getMostPathogenicVarType());
			newAnnotations.add(a);
			// hack to deal with fusion protein Jannovar nonsense - ? should the separate genes not be part of the
			// annotation anyway - don't seem to be, should maybe not do this split
			if (geneSymbol.contains("-")) {
				String[] separateGeneSymbols = geneSymbol.split("-");
				for (String separateGeneSymbol : separateGeneSymbols) {
					geneSymbols.add(separateGeneSymbol);
					variantEffects.add(VariantEffect.CUSTOM);// for - split entries do not know effect
					newAnnotations.add(null);
				}
			}
		}

		Annotation bestAnnotation = null;
		int i = 0;
		for (String geneSymbol : geneSymbols) {
			Gene gene = allGenes.get(geneSymbol);
			VariantEffect ve = variantEffects.get(i);
			if (gene != null && (gene.getPriorityResult(priorityType)) != null) {
				float geneScore = prioritiserScore(gene);
				if (geneScore > bestScore) {
					bestScore = geneScore;
					geneWithHighestPhenotypeScore = gene;
					variantEffectForTopHit = ve;
					if (newAnnotations.get(i) != null) {
						bestAnnotation = newAnnotations.get(i);
					} else {
						bestAnnotation = null;
					}
				}
			}
			i++;
		}
		// Keep original annotation if possible - used in RegFilter later on and for display
		List<Annotation> finalAnnotations = new ArrayList<>();
		if (bestAnnotation != null) {
			finalAnnotations.add(bestAnnotation);
		}

		if (prioritiserScore(currentlyAssignedGene) == bestScore) {
			// don't move the assignment if there is nowhere better to go...
			return;
		}
		if (geneWithHighestPhenotypeScore == null) {
			return;
		}
		/*
		 * Only reassign variant effect if it has not already been flagged by RegFeatureDao as a regulatory region
		 * variant. Otherwise TAD reassignment and subsequent reg feature filter fail to work as expected
		 */
		if (variantEvaluation.getVariantEffect() != VariantEffect.REGULATORY_REGION_VARIANT) {
			variantEvaluation.setVariantEffect(variantEffectForTopHit);
		}
		variantEvaluation.setAnnotations(finalAnnotations);
		variantEvaluation.setEntrezGeneId(geneWithHighestPhenotypeScore.getEntrezGeneID());
		variantEvaluation.setGeneSymbol(geneWithHighestPhenotypeScore.getGeneSymbol());
	}

}
