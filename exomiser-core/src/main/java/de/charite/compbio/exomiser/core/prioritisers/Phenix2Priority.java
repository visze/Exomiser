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

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.io.Files;

import de.charite.compbio.exomiser.core.model.Gene;
import drseb.BoqaService;
import drseb.BoqaService.ResultEntry;
import ontologizer.go.Ontology;
import ontologizer.go.Term;

/**
 * Filter variants according to the BOQA-score between given query and HPO-annotated disease genes.
 *
 * @author Sebastian Koehler
 * @version 0.01 (20 April, 2016)
 */
public class Phenix2Priority implements Prioritiser {

    private static final Logger logger = LoggerFactory.getLogger(Phenix2Priority.class);

    private static final String geneIdPrefix = "ENTREZGENE";


    /**
     * The HPO terms entered by the user describing the individual who is being
     * sequenced by exome-sequencing or clinically relevant genome panel.
     */
    private HashSet<Term> hpoQueryTerms;

	/**
	 * Proxy to BOQA
	 */
	private BoqaService boqaService;

	private int analysedGenes;

	private int numberQueryTerms;

    /**
     * Create a new instance of the Phenix2Priority.
     *
     */
    public Phenix2Priority(String hpOboFile, String geneAnnotationFile, List<String> hpoQueryTermIds) {

        /*
         * Quick sanity check of gene-annotation file:<br>
         *  - File must exist <br>
         *  - first line should start with "ENTREZGENE"<br>
         *  - first line should have 14 elements (tab-sep)
         */
        try {
			String firstLine = Files.readFirstLine(new File(geneAnnotationFile), Charset.defaultCharset());
			if (!firstLine.startsWith(geneIdPrefix+"\t")) {
				 throw new PhenixException("Lines of gene-annotation file in Phenix2 must start with "+geneIdPrefix);
			}
			if (!(firstLine.split("\t").length==14)) {
				 throw new PhenixException("Lines of gene-annotation file in Phenix2 must have 14 elements");
			}
		} catch (IOException e) {
			throw new PhenixException("I/O-Problem with gene-annotation file given to Phenix2: "+e.getMessage());
		}
        
        boqaService = new BoqaService(hpOboFile, geneAnnotationFile);
        
        if (hpoQueryTermIds.isEmpty()) {
            throw new PhenixException("Please supply some HPO terms. PhenIX2 is unable to prioritise genes without these.");
        }

        Ontology hpo = boqaService.getOntology();
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
    }

    /**
     * STUB CONSTRUCTOR - ONLY USED FOR TESTING PURPOSES TO AVOID NULL POINTERS FROM ORIGINAL CONSTRUCTOR. DO NOT USE FOR PRODUCTION CODE!!!!
     * @param hpoIds
     * @param symmetric 
     */
    public Phenix2Priority(List<String> hpoQueryTermIds) {
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

        HashMap<String, ResultEntry> scoredGenes = boqaService.scoreItems(hpoQueryTerms);
        HashMap<String, Double> geneid2score = new HashMap<>();
        for (ResultEntry resultEntry : scoredGenes.values()) {
        	geneid2score.put(resultEntry.getItemRealId(), resultEntry.getScore());
        }
        
        
        for (Gene gene : genes) {
        	String geneId = geneIdPrefix+":"+gene.getEntrezGeneID();
        	if (geneid2score.containsKey(geneId)) {
        		Phenix2PriorityResult phenix2result = new Phenix2PriorityResult(geneid2score.get(geneId));
        		gene.addPriorityResult(phenix2result);
        	}
        	else {
        		Phenix2PriorityResult phenix2result = new Phenix2PriorityResult(0);
        		gene.addPriorityResult(phenix2result);
        	}
        }
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
