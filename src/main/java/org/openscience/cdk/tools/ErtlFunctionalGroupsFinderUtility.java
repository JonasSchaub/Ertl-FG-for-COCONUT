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

import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Objects;
import java.util.stream.IntStream;

/**
 *
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class ErtlFunctionalGroupsFinderUtility {

    //<editor-fold desc="Private static final class constants">
    /**
     *
     */
    private static final int[] VALID_ATOMIC_NUMBERS = new int[] {1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86};
    //</editor-fold>

    //<editor-fold desc="Constructors">
    /**
     *
     */
    public ErtlFunctionalGroupsFinderUtility() {

    }
    //</editor-fold>

    //<editor-fold desc="Public static final properties">
    /**
     * Returns an int array containing all atomic numbers that can be passed on to ErtlFunctionalGroupsFinder.find().
     * All other atomic numbers are invalid because they represent metal, metalloid or pseudo ('R') atoms.
     *
     * @return all valid atomic numbers for ErtlFunctionalGroupsFinder.find()
     */
    public static final int[] getValidAtomicNumbers() {
        return ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS;
    }
    //</editor-fold>

    //<editor-fold desc="Public static methods">
    /**
     * Checks whether the given molecule consists of two or more unconnected structures, e.g. ion and counter-ion.
     *
     * @param aMolecule the molecule to check
     * @return true, if the molecule consists of two or more unconnected structures
     * @throws NullPointerException if the given molecule is 'null'
     */
    public static boolean isStructureUnconnected(IAtomContainer aMolecule) throws NullPointerException {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'.");
        return (!ConnectivityChecker.isConnected(aMolecule));
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
     * @throws NullPointerException if the given atom is 'null'
     */
    public static boolean isAtomCharged(IAtom anAtom) throws NullPointerException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
        Integer tmpFormalCharge = anAtom.getFormalCharge();
        return (tmpFormalCharge!= 0 && tmpFormalCharge != null);
    }

    /**
     * Checks whether a given atom is a metal, metalloid or pseudo atom judged by its atomic number. Atoms with invalid
     * atomic numbers (metal, metalloid or pseudo ('R') atoms) can not be passed on to ErtlFunctionalGroupsFinder.find()
     * but should be discarded.
     *
     * @param anAtom the atom to check
     * @return true, if the atomic number is invalid
     * @throws NullPointerException if the given atom is 'null'
     */
    public static boolean isAtomicNumberInvalid(IAtom anAtom) throws NullPointerException {
        Objects.requireNonNull(anAtom, "Given atom is 'null'.");
        Integer tmpAtomicNumber = anAtom.getAtomicNumber();
        if (tmpAtomicNumber == null) {
            return false;
        }
        boolean tmpIsAtomicNumberValid = IntStream.of(ErtlFunctionalGroupsFinderUtility.VALID_ATOMIC_NUMBERS).anyMatch(x -> x == tmpAtomicNumber);
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

}
