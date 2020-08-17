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
 * - Update to more recent version of COCONUT
 * - Find a way to print non-String value fields of the documents
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
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Collections;
import java.util.Map;
import java.util.Objects;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * TODO
 *
 * @author Jonas Schaub
 */
public class CreateCoconutSdfTest {
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
    private static final String OUTPUT_FOLDER_NAME = CreateCoconutSdfTest.class.getSimpleName() + "_Output";

    /**
     * Name of the output file
     */
    private static final String OUTPUT_File_NAME = CreateCoconutSdfTest.DATABASE_NAME + "_"
            + CreateCoconutSdfTest.COLLECTION_NAME + ".sdf";

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
    private static final String SMILES_CODE_KEY = "smiles";

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(CreateCoconutSdfTest.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Public test methods">
    @Test
    public void exportCoconutToSdf() throws Exception {
        MongoClientSettings.Builder tmpBuilder = MongoClientSettings.builder();
        ServerAddress tmpAddress = new ServerAddress(CreateCoconutSdfTest.HOST, CreateCoconutSdfTest.PORT);
        tmpBuilder.applyToClusterSettings(builder -> builder.hosts(Collections.singletonList(tmpAddress)));
        MongoClientSettings tmpSettings = tmpBuilder.build();
        MongoClient tmpMongoClient = MongoClients.create(tmpSettings);
        MongoDatabase tmpDatabase = tmpMongoClient.getDatabase(CreateCoconutSdfTest.DATABASE_NAME);
        MongoCollection<Document> tmpCollection = tmpDatabase.getCollection(CreateCoconutSdfTest.COLLECTION_NAME);
        MongoCursor<Document> tmpCursor = null;
        try {
            tmpCursor = tmpCollection.find().iterator();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            CreateCoconutSdfTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Connection to MongoDB successful.");
        System.out.println("Collection " + CreateCoconutSdfTest.COLLECTION_NAME + " in database " + CreateCoconutSdfTest.DATABASE_NAME + " is loaded.");
        String tmpOutputFolderPath = (new File(CreateCoconutSdfTest.OUTPUT_FOLDER_NAME)).getAbsolutePath() + File.separator;
        System.out.println("Output directory: " + tmpOutputFolderPath);
        File tmpOutputFile = new File(tmpOutputFolderPath + CreateCoconutSdfTest.OUTPUT_File_NAME);
        if (!tmpOutputFile.exists()) {
            tmpOutputFile.mkdirs();
        }
        FileOutputStream tmpFileOut = new FileOutputStream(tmpOutputFile);
        SDFWriter tmpWriter = new SDFWriter(tmpFileOut);
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputFolderPath + CreateCoconutSdfTest.LOG_FILE_NAME);
        } catch (IOException anIOException) {
            CreateCoconutSdfTest.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger.getLogger("").addHandler(tmpLogFileHandler);
        Logger.getLogger("").setLevel(Level.ALL);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpFilteredCounter = 0;
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
                tmpID = tmpCurrentDoc.getString(CreateCoconutSdfTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(CreateCoconutSdfTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                if (Objects.isNull(tmpMolecule)) {
                    tmpFilteredCounter++;
                    continue;
                }
                tmpMolecule.setProperties((Map) tmpCurrentDoc);
                tmpWriter.write(tmpMolecule);
                tmpCurrentDoc = null;
                tmpID = "";
                tmpSmilesCode = "";
                tmpMolecule = null;
            } catch (Exception anException) {
                CreateCoconutSdfTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Filtered counter: " + tmpFilteredCounter);
        tmpWriter.close();
        tmpLogFileHandler.flush();
        tmpCursor.close();
    }
    //</editor-fold>
}
