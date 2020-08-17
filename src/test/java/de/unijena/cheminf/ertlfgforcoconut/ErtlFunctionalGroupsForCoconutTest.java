/*
 * Uses
 * ErtlFunctionalGroupsFinder for CDK
 * to extract functional groups from a database of natural products
 * Copyright (C) 2020 Jonas Schaub
 *
 * Source code is available at <https://github.com/JonasSchaub/Ertl-FG-for-COCONUT>
 * ErtlFunctionalGroupsFinder for CDK is available at <https://github.com/zielesny/ErtlFunctionalGroupsFinder>
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
package de.unijena.cheminf.ertlfgforcoconut;

/**
 * TODO:
 * -
 */

import com.mongodb.MongoClientSettings;
import com.mongodb.MongoTimeoutException;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoClients;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.MongoDatabase;
import org.bson.Document;
import org.junit.Assume;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinder;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinderUtility;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Objects;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * This test class can be used to load a MongoDB database (or rather one collection in it) containing molecules
 * (represented as SMILES code strings), extracting all Ertl functional groups from these molecules, and compiling the
 * information which functional groups occurred in each molecule in a text (csv) file.
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class ErtlFunctionalGroupsForCoconutTest {
    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
    /**
     * Host of the MongoDB instance
     */
    private static final String HOST = "localhost";

    /**
     * Port where the MongoDB instance is running
     */
    private static final int PORT = 27017;

    /**
     * Name of the MongoDB database
     */
    private static final String DATABASE_NAME = "COCONUTfebruary20";

    /**
     * Collection from the database to load
     */
    private static final String COLLECTION_NAME = "uniqueNaturalProduct";

    /**
     * Name of the output folder
     */
    private static final String OUTPUT_FOLDER_NAME = "ErtlFunctionalGroupsForCoconut_Output";

    /**
     * Name of the output file containing the compiled information about functional groups detected in the given molecules
     */
    private static final String OUTPUT_File_NAME = "Functional_groups.txt";

    /**
     * Name of the log file for exceptions
     */
    private static final String LOG_FILE_NAME = "Log.txt";

    /**
     * Name of the document variable that contains an ID of the given molecule
     */
    private static final String ID_KEY = "coconut_id";

    /**
     * Name of the document variable that contains the SMILES code string of the given molecule
     */
    private static final String SMILES_CODE_KEY = "smiles"; //clean_smiles or smiles

    /**
     * Separator for the results file (csv)
     */
    private static final String OUTPUT_FILE_SEPARATOR = ",";

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(ErtlFunctionalGroupsForCoconutTest.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Public test methods">
    /**
     * This test method loads a MongoDB database (or rather one collection in it) containing molecules
     * (represented as SMILES code strings), extracts all Ertl functional groups from these molecules, and compiles the
     * information which functional groups occurred in each molecule in a text (csv) file.
     * <br>If no connection to MongoDB can be made, the test is ignored.
     *
     * @throws Exception if anything unexpected happens; all exceptions caused by the respective molecules are caught and logged
     */
    @Test
    public void extractFunctionalGroupsFromCoconut() throws Exception {
        MongoClientSettings.Builder tmpBuilder = MongoClientSettings.builder();
        ServerAddress tmpAddress = new ServerAddress(ErtlFunctionalGroupsForCoconutTest.HOST, ErtlFunctionalGroupsForCoconutTest.PORT);
        tmpBuilder.applyToClusterSettings(builder -> builder.hosts(Collections.singletonList(tmpAddress)));
        MongoClientSettings tmpSettings = tmpBuilder.build();
        MongoClient tmpMongoClient = MongoClients.create(tmpSettings);
        MongoDatabase tmpDatabase = tmpMongoClient.getDatabase(ErtlFunctionalGroupsForCoconutTest.DATABASE_NAME);
        MongoCollection<Document> tmpCollection = tmpDatabase.getCollection(ErtlFunctionalGroupsForCoconutTest.COLLECTION_NAME);
        MongoCursor<Document> tmpCursor = null;
        try {
            tmpCursor = tmpCollection.find().iterator();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            ErtlFunctionalGroupsForCoconutTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Connection to MongoDB successful.");
        System.out.println("Collection " + ErtlFunctionalGroupsForCoconutTest.COLLECTION_NAME + " in database " + ErtlFunctionalGroupsForCoconutTest.DATABASE_NAME + " is loaded.");
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        String tmpOutputFolderPath = (new File(tmpClassLoader.getResource(ErtlFunctionalGroupsForCoconutTest.OUTPUT_FOLDER_NAME).getFile())).getAbsolutePath() + File.separator;
        System.out.println("Output directory: " + tmpOutputFolderPath);
        File tmpOutputFile = new File(tmpOutputFolderPath + ErtlFunctionalGroupsForCoconutTest.OUTPUT_File_NAME);
        FileWriter tmpResultsWriter = null;
        try {
            tmpResultsWriter = new FileWriter(tmpOutputFile);
        } catch (IOException anIOException) {
            ErtlFunctionalGroupsForCoconutTest.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while creating the results file. Test is abandoned.");
            Assume.assumeTrue(false);
        }
        PrintWriter tmpResultsPrinter = new PrintWriter(tmpResultsWriter);
        tmpResultsPrinter.println(ErtlFunctionalGroupsForCoconutTest.ID_KEY + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR
                + "FgSMILES" + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR + "FgPseudoSMILES"
                + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        tmpResultsPrinter.flush();
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        Aromaticity tmpAromaticityModel = new Aromaticity(ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.cdkAromaticSet()));
        ErtlFunctionalGroupsFinder tmpErtlFinder = ErtlFunctionalGroupsFinderUtility.getErtlFunctionalGroupsFinderGeneralizingMode();
        MoleculeHashGenerator tmpHashGenerator = ErtlFunctionalGroupsFinderUtility.getFunctionalGroupHashGenerator();
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputFolderPath + ErtlFunctionalGroupsForCoconutTest.LOG_FILE_NAME);
        } catch (IOException anIOException) {
            ErtlFunctionalGroupsForCoconutTest.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger.getLogger("").addHandler(tmpLogFileHandler);
        Logger.getLogger("").setLevel(Level.ALL);
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        List<IAtomContainer> tmpFunctionalGroupsGeneralized;
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpFilteredCounter = 0;
        int tmpNoneDetectedCounter = 0;
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
                tmpID = tmpCurrentDoc.getString(ErtlFunctionalGroupsForCoconutTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(ErtlFunctionalGroupsForCoconutTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                if (Objects.isNull(tmpMolecule)) {
                    tmpFilteredCounter++;
                    continue;
                }
                tmpMolecule.setTitle(tmpID);
                tmpMolecule = ErtlFunctionalGroupsFinderUtility.applyFiltersAndPreprocessing(tmpMolecule, tmpAromaticityModel);
                if (Objects.isNull(tmpMolecule)) {
                    tmpResultsPrinter.println(tmpID + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR + "[got filtered]");
                    tmpResultsPrinter.flush();
                    tmpFilteredCounter++;
                    continue;
                }
                tmpFunctionalGroupsGeneralized = tmpErtlFinder.find(tmpMolecule, false);
                if (tmpFunctionalGroupsGeneralized.isEmpty()) {
                    tmpResultsPrinter.println(tmpID + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR + "[none detected]");
                    tmpResultsPrinter.flush();
                    tmpNoneDetectedCounter++;
                    continue;
                }
                HashMap<Long, IAtomContainer> tmpResultsMap = new HashMap<>(tmpFunctionalGroupsGeneralized.size(), 1);
                for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroupsGeneralized) {
                    Long tmpHashCode = tmpHashGenerator.generate(tmpFunctionalGroup);
                    if (tmpResultsMap.keySet().contains(tmpHashCode)) {
                        int tmpFrequency = tmpResultsMap.get(tmpHashCode).getProperty("FREQUENCY");
                        tmpResultsMap.get(tmpHashCode).setProperty("FREQUENCY", tmpFrequency + 1);
                    } else {
                        tmpFunctionalGroup.setProperty("FREQUENCY", 1);
                        tmpResultsMap.put(tmpHashCode, tmpFunctionalGroup);
                    }
                }
                String tmpResultsLine = tmpID;
                for (Long tmpHashCode : tmpResultsMap.keySet()) {
                    IAtomContainer tmpFunctionalGroup = tmpResultsMap.get(tmpHashCode);
                    String tmpFGSmilesCode = tmpSmiGen.create(tmpFunctionalGroup);
                    String tmpFGPseudoSmilesCode = ErtlFunctionalGroupsFinderUtility.createPseudoSmilesCode(tmpFunctionalGroup);
                    int tmpFrequency = tmpFunctionalGroup.getProperty("FREQUENCY");
                    tmpResultsLine = tmpResultsLine + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR
                            + tmpFGSmilesCode + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR
                            + tmpFGPseudoSmilesCode + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR + tmpFrequency;
                }
                tmpResultsPrinter.println(tmpResultsLine);
                tmpResultsPrinter.flush();
                tmpCurrentDoc = null;
                tmpID = "";
                tmpSmilesCode = "";
                tmpMolecule = null;
            } catch (Exception anException) {
                ErtlFunctionalGroupsForCoconutTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                tmpResultsPrinter.println(tmpID + ErtlFunctionalGroupsForCoconutTest.OUTPUT_FILE_SEPARATOR + "[exception occurred]");
                tmpResultsPrinter.flush();
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Filtered counter: " + tmpFilteredCounter);
        System.out.println("No functional groups detected: " + tmpNoneDetectedCounter);
        tmpResultsPrinter.close();
        try {
            tmpResultsWriter.close();
        } catch (IOException anIOException) {
            ErtlFunctionalGroupsForCoconutTest.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
        }
        tmpCursor.close();
    }
    //</editor-fold>
}
