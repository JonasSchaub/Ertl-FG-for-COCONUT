/**
 * Test of Utilities for
 * ErtlFunctionalGroupsFinder for CDK
 * Copyright (C) 2020 Jonas Schaub
 *
 * Source code is available at <https://github.com/JonasSchaub/Ertl-FG-for-COCONUT>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.openscience.cdk.tools;

/**
 * TODO:
 * - Implement tests for hash generator settings ands preprocessing and copy method? And other functionalities
 * - Add docs
 */

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.HashMap;

/**
 *
 */
public class ErtlFunctionalGroupsFinderUtilityTest {
    /**
     *
     * @throws Exception
     */
    @Test
    public void testPseudoSmilesGeneration() throws Exception {
        HashMap<String, String> tmpTestPairsMap = new HashMap<>(20, 1);
        tmpTestPairsMap.put("*n(*)*", "RN*(R)R");
        tmpTestPairsMap.put("*O*", "ROR");
        tmpTestPairsMap.put("[H]O[c]", "[H]O[C*]");
        tmpTestPairsMap.put("*OC(*)=O", "ROC(R)=O");
        tmpTestPairsMap.put("*o*", "RO*R");
        tmpTestPairsMap.put("[c]=O", "[C*]=O");
        tmpTestPairsMap.put("[As]", "[As]");
        tmpTestPairsMap.put("[Po]", "[Po]");
        //The CDK SmilesParser cannot parse the element Uup, it gets turned into a wildcard ('*')
        tmpTestPairsMap.put("[Uup]", "R");
        tmpTestPairsMap.put("[se]", "[Se*]");
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        tmpSmilesParser.kekulise(false);
        IAtomContainer tmpTestMolecule;
        String tmpPseudoSmilesCode;
        for (String tmpSmilesCode : tmpTestPairsMap.keySet()) {
            tmpTestMolecule = tmpSmilesParser.parseSmiles(tmpSmilesCode);
            tmpPseudoSmilesCode = ErtlFunctionalGroupsFinderUtility.getPseudoSmilesCode(tmpTestMolecule);
            System.out.println(tmpPseudoSmilesCode);
            Assert.assertEquals(tmpTestPairsMap.get(tmpSmilesCode), tmpPseudoSmilesCode);
        }
    }
}
