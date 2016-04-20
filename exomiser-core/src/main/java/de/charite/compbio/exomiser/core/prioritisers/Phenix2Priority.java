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

package de.charite.compbio.exomiser.core.prioritisers;

import de.charite.compbio.exomiser.core.model.Gene;
import de.charite.compbio.exomiser.core.prioritisers.util.ScoreDistribution;
import de.charite.compbio.exomiser.core.prioritisers.util.ScoreDistributionContainer;
import drseb.BoqaService;
import drseb.BoqaService.ResultEntry;
import hpo.HPOutils;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import ontologizer.go.OBOParser;
import ontologizer.go.OBOParserException;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermContainer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.io.Files;

import hpo.similarity.SimilarityUtilities;
import hpo.similarity.concepts.ResnikSimilarity;
import hpo.similarity.objects.InformationContentObjectSimilarity;
import sonumina.math.graph.SlimDirectedGraphView;

/**
 * Filter variants according to the BOQA-score between given query and HPO-annotated disease genes.
 *
 * @author Sebastian Koehler
 * @version 0.01 (20 April, 2016)
 */
public class Phenix2Priority implements Prioritiser {

    private static final Logger logger = LoggerFactory.getLogger(Phenix2Priority.class);

    /**
     * The HPO as Ontologizer-Ontology object
     */
    private Ontology hpo;


    /**
     * A list of error-messages
     */
    private List<String> errorMessages = null;


    /**
     * The HPO terms entered by the user describing the individual who is being
     * sequenced by exome-sequencing or clinically relevant genome panel.
     */
    private HashSet<Term> hpoQueryTerms;

    private float DEFAULT_SCORE = 0f;

    private Map<String, List<Term>> geneId2annotations;

    private final ScoreDistributionContainer scoredistributionContainer = new ScoreDistributionContainer();

    private int numberQueryTerms;
    /**
     * A counter of the number of genes that could not be found in the database
     * as being associated with a defined disease gene.
     */
    private int offTargetGenes = 0;
    /**
     * Total number of genes used for the query, including genes with no
     * associated disease.
     */
    private int analysedGenes;

    /**
     * Path to the directory that has the files needed to calculate the score
     * distribution.
     */
    private String scoredistributionFolder;
    /**
     * Keeps track of the maximum semantic similarity score to date
     */
    private double maxBOQAscore = 0d;

	private BoqaService boqaService;

    /**
     * Create a new instance of the Phenix2Priority.
     *
     */
    public Phenix2Priority(String hpOboFile, String geneAnnotationFile, List<String> hpoQueryTermIds) {

        if (hpoQueryTermIds.isEmpty()) {
            throw new PhenixException("Please supply some HPO terms. PhenIX2 is unable to prioritise genes without these.");
        }

        hpoQueryTerms = new HashSet<>();        
        for (String termIdString : hpoQueryTermIds) {
            Term t = hpo.getTermIncludingAlternatives(termIdString);
            if (t != null) {
            	hpoQueryTerms.add(t);
            } else {
                logger.error("invalid term-id given: " + termIdString);
            }
        }
        numberQueryTerms = hpoQueryTerms.size();
        
        /*
         * Quick sanity check of gene-annotation file:<br>
         *  - File must exist <br>
         *  - first line should start with "ENTREZGENE"<br>
         *  - first line should have 14 elements (tab-sep)
         */
        try {
			String firstLine = Files.readFirstLine(new File(geneAnnotationFile), Charset.defaultCharset());
			if (!firstLine.startsWith("ENTREZGENE\t")) {
				 throw new PhenixException("Lines of gene-annotation file in Phenix2 must start with ENTREZGENE");
			}
			if (!(firstLine.split("\t").length==14)) {
				 throw new PhenixException("Lines of gene-annotation file in Phenix2 must have 14 elements");
			}
		} catch (IOException e) {
			throw new PhenixException("I/O-Problem with gene-annotation file given to Phenix2: "+e.getMessage());
		}
        
        boqaService = new BoqaService(hpOboFile, geneAnnotationFile);
    }

    /**
     * Flag to output results of filtering with BOQA.
     */
    @Override
    public PriorityType getPriorityType() {
        return PriorityType.PHENIX2_PRIORITY;
    }

    /**
     * Prioritize a list of candidate {@link exomizer.exome.Gene Gene} objects
     * (the candidate genes have rare, potentially pathogenic variants).
     *
     * @param genes List of candidate genes.
     * @see exomizer.filter.Filter#filter_list_of_variants(java.util.ArrayList)
     */
    @Override
    public void prioritizeGenes(List<Gene> genes) {
        analysedGenes = genes.size();

        // TODO ... run boqa and set scores for each Gene-Object
        HashMap<String, ResultEntry> scoredGenes = boqaService.scoreItems(hpoQueryTerms);
        for (Gene gene : genes) {
            Phenix2PriorityResult phenomizerRelScore = scoreVariantHPO(gene);
            gene.addPriorityResult(phenomizerRelScore);
        }
        normalizeBoqaScores(genes);
    }

    /**
     * The gene relevance scores are to be normalized to lie between zero and
     * one. This function, which relies upon the variable {@link #maxBOQAscore}
     * being set in {@link #scoreVariantHPO(Gene)}, divides each score by
     * {@link #maxBOQAscore}, which has the effect of putting the BOQA scores
     * in the range [0..1].
     */
    private void normalizeBoqaScores(List<Gene> genes) {
        if (maxBOQAscore < 1) {
            return;
        }
        Phenix2PriorityResult.setNormalizationFactor(1d / maxBOQAscore);
    }

    /**
     * @param gene A {@link exomizer.exome.Gene Gene} whose score is to be
     * determined.
     */
    private Phenix2PriorityResult scoreVariantHPO(Gene gene) {

    	System.err.println("to be implemented!");
    	return null;
//        int entrezGeneId = gene.getEntrezGeneID();
//        String entrezGeneIdString = entrezGeneId + "";
//
//        if (!geneId2annotations.containsKey(entrezGeneIdString)) {
//            //System.err.println("INVALID GENE GIVEN (will set to default-score): Entrez ID: " + g.getEntrezGeneID() + " / " + g.getGeneSymbol());
//            this.offTargetGenes++;
//            return new PhenixPriorityResult(DEFAULT_SCORE);
//        }
//
//        List<Term> annotationsOfGene = geneId2annotations.get(entrezGeneIdString);
//
//        double similarityScore = similarityMeasure.computeObjectSimilarity( (ArrayList<Term>) hpoQueryTerms, (ArrayList<Term>) annotationsOfGene);
//        if (similarityScore > maxSemSim) {
//            maxSemSim = similarityScore;
//        }
//        if (Double.isNaN(similarityScore)) {
//            errorMessages.add("Error: score was NAN for gene:" + gene + " : " + hpoQueryTerms + " <-> " + annotationsOfGene);
//        }
//
//        ScoreDistribution scoreDist = scoredistributionContainer.getDistribution(entrezGeneIdString, numberQueryTerms, symmetric,
//                scoredistributionFolder);
//
//	// get the pvalue
//        double rawPvalue;
//        if (scoreDist == null) {
//            return new PhenixPriorityResult(DEFAULT_SCORE);
//        } else {
//            rawPvalue = scoreDist.getPvalue(similarityScore, 1000.);
//            rawPvalue = Math.log(rawPvalue) * -1.0; /* Negative log of p value : most significant get highest score */
//
//            if (rawPvalue > maxNegLogP) {
//                maxNegLogP = rawPvalue;
//            }
//        }
//
//        return new PhenixPriorityResult(rawPvalue, similarityScore);
    }


    private static class PhenixException extends RuntimeException {

        private PhenixException(String message) {
            super(message);
        }
    }


    @Override
    public int hashCode() {
        int hash = 7;
        hash = 13 * hash + Objects.hashCode(this.hpoQueryTerms);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Phenix2Priority other = (Phenix2Priority) obj;
        if (!Objects.equals(this.hpoQueryTerms, other.hpoQueryTerms)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return "Phenix2Priority{" + "hpoQueryTerms=" + hpoQueryTerms + '}';
    }

}
