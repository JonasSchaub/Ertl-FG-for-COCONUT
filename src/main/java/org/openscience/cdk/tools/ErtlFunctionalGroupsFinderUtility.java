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

import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

import java.util.Objects;
import java.util.stream.IntStream;

/**
 * This class gives utility methods for using ErtlFunctionalGroupsFinder.
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class ErtlFunctionalGroupsFinderUtility {
    //<editor-fold desc="Private static final class constants">
    /**
     * Atomic numbers that ErtlFunctionalGroupsFinder accepts, see getValidAtomicNumbers()
     */
    private static final int[] VALID_ATOMIC_NUMBERS = new int[] {1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86};
    //</editor-fold>
    //
    //<editor-fold desc="Constructors">
    //TODO: Implement non-static methods and set up logging files in this constructor?
    /**
     * Constructor (currently empty)
     */
    public ErtlFunctionalGroupsFinderUtility() {

    }
    //</editor-fold>
    //
    //<editor-fold desc="Public static final properties">
    /**
     * Returns an integer array containing all atomic numbers that can be passed on to ErtlFunctionalGroupsFinder.find().
     * All other atomic numbers are invalid because they represent metal, metalloid or pseudo ('R') atoms.
     *
     * @return all valid atomic numbers for ErtlFunctionalGroupsFinder.find()
     */
    public static final int[] getValidAtomicNumbers() {
        return ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS;
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
     * method would still accept these molecules but it is not recommended to pass them on.
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
        boolean tmpIsAtomicNumberValid = IntStream.of(ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS).anyMatch(x -> x == tmpAtomicNumberInt);
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

    //TODO: Add doc!
    /**
     *
     */
    public static IAtom neutralizeCharges(IAtom anAtom, IAtomContainer aParentMolecule) throws NullPointerException, CDKException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
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
    //</editor-fold>
    //</editor-fold>
}
