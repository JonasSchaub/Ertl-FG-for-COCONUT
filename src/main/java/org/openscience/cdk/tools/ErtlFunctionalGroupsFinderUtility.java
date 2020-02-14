/**
 * Utilities for
 * ErtlFunctionalGroupsFinder for CDK
 * Copyright (C) 2020 Jonas Schaub
 *
 * Source code is available at <https://github.com/zielesny/ErtlFunctionalGroupsFinder>
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
 * - Implement test class for hash generator settings ands preprocessing and copy method?
 * - Add check for valid atomic number in every method(?)
 * - add note in docs that given params will be changed and copy() should be used if that is not desired (also note that
 * the exact same objects are returned)
 * - add method for FG frequency calculation requiring an iterator?
 * - Add logging methods using the class logger? Log molecule info alongside already logged exceptions?
 */

import org.openscience.cdk.Atom;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.hash.AtomEncoder;
import org.openscience.cdk.hash.BasicAtomEncoder;
import org.openscience.cdk.hash.HashGeneratorMaker;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This class gives utility methods for using ErtlFunctionalGroupsFinder.
 * <br>NOTE: It is not implemented having parallelized operations in mind!
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public final class ErtlFunctionalGroupsFinderUtility {
    //<editor-fold desc="Private static final class constants">
    //<editor-fold desc="General">
    /**
     * Atomic numbers that ErtlFunctionalGroupsFinder accepts, see getValidAtomicNumbers()
     */
    private static final int[] VALID_ATOMIC_NUMBERS = new int[] {1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86};

    /**
     * Atomic numbers that ErtlFunctionalGroupsFinder accepts, loaded into a hash set for quick determination; set is
     * filled in static initializer (see below)
     */
    private static final HashSet<Integer> VALID_ATOMIC_NUMBERS_SET = new HashSet<>(20, 1);

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(ErtlFunctionalGroupsFinderUtility.class.getName());
    //</editor-fold>
    //
    //<editor-fold defaultstate="collapsed" desc="Pseudo SMILES code">
    /**
     * Type of the generated SMILES codes
     */
    private static final int SMILES_GENERATOR_OUTPUT_MODE = SmiFlavor.Unique;

    /**
     * SmilesGenerator for generating SMILES and pseudo SMILES representations
     */
    private static final SmilesGenerator SMILES_GENERATOR = new SmilesGenerator(ErtlFunctionalGroupsFinderUtility.SMILES_GENERATOR_OUTPUT_MODE);

    /**
     * A map that gives a certain element symbol for a placeholder atom marking a specific aromatic atom in pseudo SMILES
     * creation
     */
    private static final HashMap<String, String> PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP = new HashMap<>(10, 1);

    /**
     * A map that gives the pseudo SMILES representation for a specific placeholder element from
     * pseudoSmilesAromaticElementToPlaceholderElementMap
     */
    private static final HashMap<String, String> PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP = new HashMap<>(10, 1);

    /**
     * Pseudo SMILES representation of an aromatic C atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_CARBON = "C*";

    /**
     * Pseudo SMILES representation of an aromatic N atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_NITROGEN = "N*";

    /**
     * Pseudo SMILES representation of an aromatic S atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_SULPHUR = "S*";

    /**
     * Pseudo SMILES representation of an aromatic O atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_OXYGEN = "O*";

    /**
     * Pseudo SMILES representation of an aromatic Se atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_SELENIUM = "Se*";

    /**
     * Pseudo SMILES representation of an aromatic P atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_PHOSPHOR = "P*";

    /**
     * Pseudo SMILES representation of an undefined pseudo atom
     */
    private static final String PSEUDO_SMILES_R_ATOM = "R";
    //</editor-fold>
    //</editor-fold>
    //
    //<editor-fold desc="Static initializer">
    /**
     * TODO: Add doc
     */
    static {
        for (int i : ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS) {
            ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS_SET.add(i);
        }
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("C", "Ce");
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("N", "Nd");
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("S", "Sm");
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("O", "Os");
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("Se", "Sc");
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("P", "Pm");
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.put("R", "Es");

        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Es",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_R_ATOM);
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Pm",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_PHOSPHOR);
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Sc",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_SELENIUM);
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Os",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_OXYGEN);
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Sm",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_SULPHUR);
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Nd",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_NITROGEN);
        ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.put("Ce",
                ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_CARBON);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    //TODO: Implement non-static methods and set up logging files in this constructor? Or static initializer? Discard constructor?
    /**
     * Constructor (currently empty)
     */
    public ErtlFunctionalGroupsFinderUtility() {

    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static properties">
    /**
     * Returns an integer array containing all atomic numbers that can be passed on to ErtlFunctionalGroupsFinder.find().
     * All other atomic numbers are invalid because they represent metal, metalloid or pseudo ('R') atoms.
     *
     * @return all valid atomic numbers for ErtlFunctionalGroupsFinder.find()
     */
    public static int[] getValidAtomicNumbers() {
        return Arrays.copyOf(ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS,
                ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS.length);
    }

    /**
     * TODO: Add doc
     */
    public static MoleculeHashGenerator getFunctionalGroupHashGenerator() {
        MoleculeHashGenerator tmpHashGenerator = new HashGeneratorMaker()
                .depth(8)
                .elemental()
                /*following line is used instead of .orbital() because the atom hybridizations take more information into
                account than the bond order sum but that is not required here*/
                /*Note: This works here because the ErtlFunctionalGroupsFinder extracts the relevant atoms and bonds only
                resulting in incomplete valences that can be used here in this way*/
                .encode(BasicAtomEncoder.BOND_ORDER_SUM)
                .encode(CustomAtomEncoder.AROMATICITY) //See enum CustomAtomEncoder below
                .molecular();
        return tmpHashGenerator;
    }

    /**
     * TODO: Add doc
     */
    public static ErtlFunctionalGroupsFinder getErtlFunctionalGroupsFinderGeneralizingMode() {
        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.DEFAULT);
        return tmpEFGF;
    }

    /**
     * TODO: Add doc
     */
    public static ErtlFunctionalGroupsFinder getErtlFunctionalGroupsFinderNonGeneralizing() {
        ErtlFunctionalGroupsFinder tmpEFGF = new ErtlFunctionalGroupsFinder(ErtlFunctionalGroupsFinder.Mode.NO_GENERALIZATION);
        return tmpEFGF;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static methods">
    //<editor-fold desc="Queries for filtering">
    /**
     * Checks whether the given molecule consists of two or more unconnected structures, e.g. ion and counter-ion.
     *
     * @param aMolecule the molecule to check
     * @return true, if the molecule consists of two or more unconnected structures
     * @throws NullPointerException if the given molecule is 'null'
     */
    public static boolean isStructureUnconnected(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        boolean tmpIsConnected = ConnectivityChecker.isConnected(aMolecule);
        return (!tmpIsConnected);
    }

    /**
     * Checks whether the atom count or bond count of the given molecule is zero. The ErtlFunctionalGroupsFinder.find()
     * method would still accept these molecules but it is not recommended to pass them on (simply makes not much sense).
     *
     * @param aMolecule the molecule to check
     * @return true, if the atom or bond count of the molecule is zero
     * @throws NullPointerException if the given molecule is 'null'
     */
    public static boolean isAtomOrBondCountZero(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int tmpAtomCount = aMolecule.getAtomCount();
        int tmpBondCount = aMolecule.getBondCount();
        return (tmpAtomCount == 0 || tmpBondCount == 0);
    }

    /**
     * Iterates through all atoms in the given molecule and checks whether they are charged. If this method returns
     * 'true', the molecule can not be passed on to ErtlFunctionalGroupsFinder.find() but should be discarded or the
     * charges neutralized.<br>
     *     If no charged atoms are found, this method scales linearly with O(n) with n: number of atoms in the given
     *     molecule.
     *
     * @param aMolecule the molecule to check
     * @return true, if the molecule contains one or more charged atoms
     * @throws NullPointerException if the given molecule (or one of its atoms) is 'null'
     */
    public static boolean isMoleculeCharged(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int tmpAtomCount = aMolecule.getAtomCount();
        if (tmpAtomCount == 0) {
            return false;
        }
        Iterable<IAtom> tmpAtoms = aMolecule.atoms();
        boolean tmpIsAtomCharged;
        for (IAtom tmpAtom : tmpAtoms) {
            //Throws NullPointerException if tmpAtom is 'null'
            tmpIsAtomCharged = ErtlFunctionalGroupsFinderUtility.isAtomCharged(tmpAtom);
            if (tmpIsAtomCharged) {
                return true;
            }
        }
        return false;
    }

    /**
     * Checks whether a given atom is charged.
     *
     * @param anAtom the atom to check
     * @return true, if the atom is charged
     * @throws NullPointerException if the given atom is 'null' or its formal charge
     */
    public static boolean isAtomCharged(IAtom anAtom) throws NullPointerException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
        Integer tmpFormalCharge = anAtom.getFormalCharge();
        Objects.requireNonNull(tmpFormalCharge, "Formal charge is 'null'.");
        return (tmpFormalCharge.intValue() != 0);
    }

    /**
     * Checks whether a given atom is a metal, metalloid or pseudo atom judged by its atomic number. Atoms with invalid
     * atomic numbers (metal, metalloid or pseudo ('R') atoms) can not be passed on to ErtlFunctionalGroupsFinder.find()
     * but should be discarded.
     *
     * @param anAtom the atom to check
     * @return true, if the atomic number is invalid or 'null'
     * @throws NullPointerException if the given atom is 'null' or its atomic number
     */
    public static boolean isAtomicNumberInvalid(IAtom anAtom) throws NullPointerException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
        Integer tmpAtomicNumber = anAtom.getAtomicNumber();
        Objects.requireNonNull(tmpAtomicNumber, "Atomic number is 'null'.");
        int tmpAtomicNumberInt = tmpAtomicNumber.intValue();
        //boolean tmpIsAtomicNumberValid = IntStream.of(ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS).anyMatch(x -> x == tmpAtomicNumberInt);
        boolean tmpIsAtomicNumberValid = ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS_SET.contains(tmpAtomicNumberInt);
        return !tmpIsAtomicNumberValid;
    }

    /**
     * Iterates through all atoms in the given molecule and checks whether their atomic numbers are invalid. If this
     * method returns 'true', the molecule can not be passed on to ErtlFunctionalGroupsFinder.find() but should be
     * discarded.<br>
     *     If no invalid atoms are found, this method scales linearly with O(n) with n: number of atoms in the given
     *     molecule.
     *
     * @param aMolecule the molecule to check
     * @return true, if the molecule contains one or more atoms with invalid atomic numbers
     * @throws NullPointerException if the given molecule (or one of its atoms) is 'null'
     */
    public static boolean containsInvalidAtomicNumbers(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        int tmpAtomCount = aMolecule.getAtomCount();
        if (tmpAtomCount == 0) {
            return false;
        }
        Iterable<IAtom> tmpAtoms = aMolecule.atoms();
        boolean tmpIsAtomicNumberInvalid;
        for (IAtom tmpAtom : tmpAtoms) {
            //Throws NullPointerException if tmpAtom is 'null'
            tmpIsAtomicNumberInvalid = ErtlFunctionalGroupsFinderUtility.isAtomicNumberInvalid(tmpAtom);
            if (tmpIsAtomicNumberInvalid) {
                return true;
            }
        }
        return false;
    }

    /**
     * TODO: Add doc
     */
    public static boolean shouldBeFiltered(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpShouldBeFiltered;
        try {
            tmpShouldBeFiltered = (ErtlFunctionalGroupsFinderUtility.containsInvalidAtomicNumbers(aMolecule)
                    || ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(aMolecule));
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.WARNING, anException.toString(), anException);
            tmpShouldBeFiltered = true;
        }
        return tmpShouldBeFiltered;
    }

    /**
     * TODO: Add doc
     */
    public static boolean shouldBePreprocessed(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpNeedsPreprocessing;
        try {
            tmpNeedsPreprocessing = (ErtlFunctionalGroupsFinderUtility.isMoleculeCharged(aMolecule)
                    || ErtlFunctionalGroupsFinderUtility.isStructureUnconnected(aMolecule));
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.WARNING, anException.toString(), anException);
            throw new NullPointerException("An unknown error occurred.");
        }
        return tmpNeedsPreprocessing;
    }

    /**
     * TODO: Add doc
     */
    public static boolean isValidArgumentForFindMethod(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is null.");
        boolean tmpIsValid;
        try {
            tmpIsValid = !(ErtlFunctionalGroupsFinderUtility.containsInvalidAtomicNumbers(aMolecule)
                    || ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(aMolecule)
                    || ErtlFunctionalGroupsFinderUtility.isMoleculeCharged(aMolecule)
                    || ErtlFunctionalGroupsFinderUtility.isStructureUnconnected(aMolecule));
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            tmpIsValid = false;
        }
        return tmpIsValid;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Preprocessing methods">
    /**
     * Returns the biggest unconnected component/structure of the given atom container, judging by the atom count. To
     * pre-check whether the atom container consists of multiple unconnected components, use isStructureUnconnected().
     * All set properties of aMolecule will be set as properties of the returned atom container.<br>
     *     Iterates through all unconnected components in the given atom container, so the method scales linearly with
     *     O(n) with n: number of unconnected components.
     *
     * @param aMolecule the molecule whose biggest unconnected component should be found
     * @return the biggest (judging by the atom count) unconnected component of the given atom container
     * @throws NullPointerException if aMolecule is null or the biggest component
     */
    public static IAtomContainer selectBiggestUnconnectedComponent(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecules is 'null'.");
        IAtomContainerSet tmpUnconnectedComponentsSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpBiggestComponent = null;
        for (IAtomContainer tmpComponent : tmpUnconnectedComponentsSet.atomContainers()) {
            if (Objects.isNull(tmpBiggestComponent) || tmpBiggestComponent.getAtomCount() < tmpComponent.getAtomCount()) {
                tmpBiggestComponent = tmpComponent;
            }
        }
        Objects.requireNonNull(tmpBiggestComponent, "The resulting biggest component is 'null'.");
        tmpBiggestComponent.setProperties(aMolecule.getProperties());
        return tmpBiggestComponent;
    }

    /**
     * Neutralizes charged atoms in the given atom container by zeroing the formal atomic charges and filling up free
     * valences with implicit hydrogen atoms (according to the CDK atom types). This procedure allows a more general
     * charge treatment than a pre-defined transformation list but may produce “wrong” structures, e.g. it turns a
     * nitro NO2 group into a structure represented by the SMILES code “[H]O[N](=O)*” with an uncharged four-bonded
     * nitrogen atom (other examples are “*[N](*)(*)*”, “[C]#[N]*” or “*S(*)(*)*”). Thus an improved charge
     * neutralization scheme is desirable for future implementations. <br>
     *     Iterates through all atoms in the given atom container, so the method scales linearly with
     *     O(n) with n: number of atoms.
     *
     * @param aMolecule the molecule to be neutralized
     * @return the same IAtomContainer instance as aMolecule but with neutralized charges
     * @throws NullPointerException if aMolecule is 'null' or one of its atoms
     * @throws CDKException if no matching atom type can be determined for one atom or there is a problem with adding
     * the implicit hydrogen atoms.
     */
    public static IAtomContainer neutralizeCharges(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        Iterable<IAtom> tmpAtoms = aMolecule.atoms();
        for (IAtom tmpAtom : tmpAtoms) {
                tmpAtom = ErtlFunctionalGroupsFinderUtility.neutralizeCharges(tmpAtom, aMolecule);
        }
        return aMolecule;
    }

    /**
     * Neutralizes a charged atom in the given parent atom container by zeroing the formal atomic charge and filling up free
     * valences with implicit hydrogen atoms (according to the CDK atom types).
     *
     * @param anAtom the atom to be neutralized
     * @param aParentMolecule the molecule the atom belongs to
     * @return the same IAtom instance as anAtom but with neutralized charges
     * @throws NullPointerException if anAtom or aParentMolecule is 'null'
     * @throws CDKException if the atom is not part of the molecule or no matching atom type can be determined for the
     * atom or there is a problem with adding the implicit hydrogen atoms.
     * @see #neutralizeCharges(IAtomContainer)
     */
    public static IAtom neutralizeCharges(IAtom anAtom, IAtomContainer aParentMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
        Objects.requireNonNull(aParentMolecule, "Given parent molecule is 'null'.");
        boolean tmpIsAtomInMolecule = aParentMolecule.contains(anAtom);
        if (!tmpIsAtomInMolecule) {
            throw new CDKException("Given atom is not part of the given atom container.");
        }
        IAtom tmpAtom = anAtom;
        Integer tmpFormalChargeObject = tmpAtom.getFormalCharge();
        if (Objects.isNull(tmpFormalChargeObject)) {
            return tmpAtom;
        }
        int tmpFormalCharge = tmpFormalChargeObject.intValue();
        if (tmpFormalCharge != 0) {
            tmpAtom.setFormalCharge(0);
            IChemObjectBuilder tmpBuilder = aParentMolecule.getBuilder();
            if (Objects.isNull(tmpBuilder)) {
                throw new CDKException("Builder of the given atom container is 'null'.");
            }
            CDKHydrogenAdder tmpHAdder = CDKHydrogenAdder.getInstance(tmpBuilder);
            CDKAtomTypeMatcher tmpMatcher = CDKAtomTypeMatcher.getInstance(tmpBuilder);
            //Can throw CDKException
            IAtomType tmpMatchedType = tmpMatcher.findMatchingAtomType(aParentMolecule, tmpAtom);
            if (Objects.isNull(tmpMatchedType)) {
                throw new CDKException("Matched atom type is 'null'.");
            }
            AtomTypeManipulator.configure(tmpAtom, tmpMatchedType);
            //Can throw CDKException
            tmpHAdder.addImplicitHydrogens(aParentMolecule, tmpAtom);
        }
        return tmpAtom;
    }

    /**
     * TODO: Add doc
     */
    public static IAtomContainer applyFiltersAndPreprocessing(IAtomContainer aMolecule, Aromaticity anAromaticityModel) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given atom container is 'null'.");
        Objects.requireNonNull(anAromaticityModel, "Given aromaticity model is 'null'.");
        IAtomContainer tmpMolecule = aMolecule;
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            //Filter
            boolean tmpIsAtomOrBondCountZero = ErtlFunctionalGroupsFinderUtility.isAtomOrBondCountZero(tmpMolecule);
            if (tmpIsAtomOrBondCountZero) {
                return null;
            }
            //From structures containing two or more unconnected structures (e.g. ions) choose the largest structure
            boolean tmpIsUnconnected = ErtlFunctionalGroupsFinderUtility.isStructureUnconnected(tmpMolecule);
            if (tmpIsUnconnected) {
                tmpMolecule = ErtlFunctionalGroupsFinderUtility.selectBiggestUnconnectedComponent(tmpMolecule);
            }
            //Filter
            boolean tmpContainsInvalidAtoms = ErtlFunctionalGroupsFinderUtility.containsInvalidAtomicNumbers(tmpMolecule);
            if (tmpContainsInvalidAtoms) {
                return null;
            }
            //Neutralize charges if there are any
            boolean tmpIsCharged = ErtlFunctionalGroupsFinderUtility.isMoleculeCharged(tmpMolecule);
            if (tmpIsCharged) {
                tmpMolecule = ErtlFunctionalGroupsFinderUtility.neutralizeCharges(tmpMolecule);
            }
            //Application of aromaticity model
            anAromaticityModel.apply(tmpMolecule);
        } catch (Exception anException) {
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.SEVERE, anException.toString(), anException);
            return null;
        }
        return tmpMolecule;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Other">
    /**
     * Aims at doing a deep copy of the given atom container, i.e. all information stored in the object is copied exactly
     * but original and copy do not share any references. The method used here writes an SD representation (via CDK's
     * SDFWriter) of the given atom container as a string and then constructs a new atom container by reading this string
     * using IteratingSDFReader. Furthermore, all properties that were stored in the original atom container are transferred
     * to the clone.
     * <p>
     *     RESTRICTIONS:
     *     - All chemical information that an SDF can represent is copied but not all 'object information', e.g. the new
     *       IAtomContainer object and its internal objects like IAtom or IBond objects will have different hash codes
     *       because these are calculated based on their memory address.
     *     - Properties stored on internal objects like IAtom or IBond objects will not be copied
     *     - Because the atom container's properties are simply transferred, copy and original will still share references
     *       to objects used as description or property
     *     - In the instantiation of the new IAtomContainer object things like unsaturated bonds might be lost
     * </p>
     *
     * @param aMolecule the atom container to copy
     * @return a deep copy of the given atom container (see restrictions)
     * @throws NullPointerException if the given atom container is 'null'
     * @throws CDKException if the SDFWriter cannot write the given atom container
     */
    public static IAtomContainer copy(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        StringWriter tmpStringWriter = new StringWriter();
        SDFWriter tmpSDFWriter = new SDFWriter(tmpStringWriter);
        boolean tmpAcceptsClass = tmpSDFWriter.accepts(aMolecule.getClass());
        if (!tmpAcceptsClass) {
            throw new CDKException("Given IChemObject/IAtomContainer implementing class can not be copied");
        }
        //Might throw CDKException
        tmpSDFWriter.write(aMolecule);
        tmpStringWriter.flush();
        try {
            tmpStringWriter.close();
            tmpSDFWriter.close();
        } catch (IOException anIOException) {
            //Should not happen because nothing is written to file
            ErtlFunctionalGroupsFinderUtility.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
        }
        String tmpSDFRepresentation = tmpStringWriter.toString();
        StringReader tmpStringReader = new StringReader(tmpSDFRepresentation);
        IteratingSDFReader tmpSDFReader = new IteratingSDFReader(tmpStringReader, aMolecule.getBuilder(), true);
        IAtomContainer tmpCopy = tmpSDFReader.next();
        //SD representation should contain string-type properties; so for the non-string-type properties, the following is done:
        Map<Object, Object> tmpProperties = aMolecule.getProperties();
        for (Object tmpKey : tmpProperties.keySet()) {
            boolean tmpCloneHasProperty = !Objects.isNull(tmpCopy.getProperty(tmpKey));
            if (!tmpCloneHasProperty) {
                tmpCopy.setProperty(tmpKey, tmpProperties.get(tmpKey));
            }
        }
        return tmpCopy;
    }

    /**
     * TODO: Add doc
     */
    public static String getPseudoSmilesCode(IAtomContainer aMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        IAtomContainer tmpMolecule = aMolecule;
        Iterable<IAtom> tmpAtoms = tmpMolecule.atoms();
        HashMap<IAtom, IAtom> tmpMapForResubstitution = new HashMap(20, 0.8f);
        for (IAtom tmpAtom: tmpAtoms) {
            boolean tmpIsAromatic = tmpAtom.isAromatic();
            boolean tmpIsPseudoAtom = (tmpAtom instanceof IPseudoAtom && "R".equals(((IPseudoAtom)tmpAtom).getLabel()));
            if (tmpIsAromatic && !tmpIsPseudoAtom) {
                String tmpSymbol = tmpAtom.getSymbol();
                if (Objects.isNull(tmpSymbol)) {
                    continue;
                }
                boolean tmpContainsKey = ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.containsKey(tmpSymbol);
                if (tmpContainsKey) {
                    String tmpReplacementElementSymbol = ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.get(tmpSymbol);
                    IAtom tmpReplacementAtom = new Atom(tmpReplacementElementSymbol);
                    Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                    //TODO: Get returned boolean and throw exception if replacement could not be made?
                    AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpAtom, tmpReplacementAtom);
                    tmpReplacementAtom.setImplicitHydrogenCount(tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
                    tmpMapForResubstitution.put(tmpReplacementAtom, tmpAtom);
                }
            }
            if (tmpIsPseudoAtom) {
                String tmpReplacementElementSymbol = ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_AROMATIC_ELEMENT_TO_PLACEHOLDER_ELEMENT_MAP.get("R");
                IAtom tmpReplacementAtom = new Atom(tmpReplacementElementSymbol);
                Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                //TODO: Get returned boolean and throw exception if replacement could not be made?
                AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpAtom, tmpReplacementAtom);
                tmpReplacementAtom.setImplicitHydrogenCount(tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
                tmpMapForResubstitution.put(tmpReplacementAtom, tmpAtom);
            }
        }
        //Might throw CDKException
        String tmpPseudoSmilesCode = ErtlFunctionalGroupsFinderUtility.SMILES_GENERATOR.create(tmpMolecule);
        for (String tmpPlaceholderElementSymbol : ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.keySet()) {
            tmpPseudoSmilesCode = tmpPseudoSmilesCode.replaceAll("(\\[" + tmpPlaceholderElementSymbol + "\\])",
                    ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.get(tmpPlaceholderElementSymbol))
                    .replaceAll("(" + tmpPlaceholderElementSymbol + ")",
                            ErtlFunctionalGroupsFinderUtility.PSEUDO_SMILES_PLACEHOLDER_ELEMENT_TO_PSEUDO_SMILES_SYMBOL_MAP.get(tmpPlaceholderElementSymbol));
        }
        for (IAtom tmpReplacementAtom: tmpMapForResubstitution.keySet()) {
            //TODO: Get returned boolean and throw exception if replacement could not be made?
            IAtom tmpOriginalAtom = tmpMapForResubstitution.get(tmpReplacementAtom);
            Integer tmpImplicitHydrogenCount = tmpReplacementAtom.getImplicitHydrogenCount();
            AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpReplacementAtom, tmpOriginalAtom);
            tmpOriginalAtom.setImplicitHydrogenCount(tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
            tmpMapForResubstitution.remove(tmpReplacementAtom, tmpOriginalAtom);
        }
        return tmpPseudoSmilesCode;
    }
    //</editor-fold>
    //</editor-fold>
}

//<editor-fold defaultstate="collapsed" desc="Enum CustomAtomEncoder">
/**
 * Custom Enumeration of atom encoders for seeding atomic hash codes.
 *
 * @author Jonas Schaub
 * @see BasicAtomEncoder
 * @see AtomEncoder
 */
enum CustomAtomEncoder implements AtomEncoder {

    /**
     * Encode whether an atom is aromatic or not. This specification is necessary to distinguish functional groups with
     * aromatic environments and those without. For example: [H]O[C] and [H]OC* (pseudo SMILES codes) should be
     * assigned different hash codes by the MoleculeHashGenerator.
     *
     * @see IAtom#isAromatic()
     */
    AROMATICITY {
        /**
         *{@inheritDoc}
         */
        @Override
        public int encode(IAtom anAtom, IAtomContainer aContainer) {
            return anAtom.isAromatic()? 3 : 2;
        }
    };
}
//</editor-fold>
