/*
 * The Exomiser - A tool to annotate and prioritize variants
 *
 * Copyright (C) 2012 - 2016  Charite Universitätsmedizin Berlin and Genome Research Ltd.
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
package de.charite.compbio.exomiser.core.dao;

import de.charite.compbio.exomiser.core.model.Variant;
import de.charite.compbio.exomiser.core.model.pathogenicity.RemmScore;
import de.charite.compbio.exomiser.core.model.pathogenicity.PathogenicityData;
import de.charite.compbio.jannovar.annotation.VariantEffect;
import htsjdk.tribble.readers.TabixReader;

import java.io.IOException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.cache.annotation.Cacheable;
import org.springframework.stereotype.Component;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
@Component
public class ESEDao {

    private final Logger logger = LoggerFactory.getLogger(ESEDao.class);

    @Autowired
    public ESEDao() {
    }

    @Cacheable(value = "ese", key = "#variant.chromosomalVariant")
    public PathogenicityData getPathogenicityData(Variant variant) {
        // ESE can only be used on Missense and Synonymous
        if (variant.getVariantEffect() != VariantEffect.MISSENSE_VARIANT || variant.getVariantEffect() != VariantEffect.SYNONYMOUS_VARIANT) {
            return new PathogenicityData();
        }
        return processResults(variant);
    }

    private PathogenicityData processResults(Variant variant) {
        String chromosome = variant.getChromosomeName();
        int start = variant.getPosition();
        int end = calculateEndPosition(variant);
        return getRemmData(chromosome, start, end);
    }

    private int calculateEndPosition(Variant variant) {
        int end = variant.getPosition();
        //these end positions are calculated according to recommendation by Max and Peter who produced the REMM score
        //don't change this unless they say. 
        if (isDeletion(variant)) {
            // test all deleted bases
            end += variant.getRef().length();
        } else if (isInsertion(variant)) {
            // test bases either side of insertion
            end += 1;
        }
        return end;
    }

    private static boolean isDeletion(Variant variant) {
        return variant.getAlt().equals("-");
    }

    private static boolean isInsertion(Variant variant) {
        return variant.getRef().equals("-");
    }
    
    private PathogenicityData getRemmData(String chromosome, int start, int end) throws NumberFormatException {
        try {
            float ncds = Float.NaN;
            String line;
//            logger.info("Running tabix with " + chromosome + ":" + start + "-" + end);
            TabixReader.Iterator results = remmTabixReader.query(chromosome + ":" + start + "-" + end);
            while ((line = results.next()) != null) {
                String[] elements = line.split("\t");
                if (Float.isNaN(ncds)) {
                    ncds = Float.parseFloat(elements[2]);
                } else {
                    ncds = Math.max(ncds, Float.parseFloat(elements[2]));
                }
            }
            //logger.info("Final score " + ncds);
            if (!Float.isNaN(ncds)) {
                return new PathogenicityData(new RemmScore(ncds));
            }
        } catch (IOException e) {
            logger.error("Unable to read from REMM tabix file {}", remmTabixReader.getSource(), e);
        }
        return new PathogenicityData();
    }

}
