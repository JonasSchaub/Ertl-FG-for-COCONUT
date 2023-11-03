# Ertl Functional Groups for the COlleCtion of Open NatUral producTs (COCONUT)
## Ertl-FG-for-COCONUT
:warning:Legacy code, not actievely maintained:warning:
<br>The algorithm for automated functional groups detection and extraction of organic molecules 
[developed by Peter Ertl](https://doi.org/10.1186/s13321-017-0225-z) and implemented on the basis of the 
Chemistry Development Kit (CDK) by Fritsch et al. ([ErtlFunctionalGroupsFinder](https://doi.org/10.1186/s13321-019-0361-8)) 
is used to supply information about functional groups for the COlleCtion of Open NatUral producTs (COCONUT) [compiled 
by Sorokina and Steinbeck](https://www.preprints.org/manuscript/201912.0332/v1).

The package <i>org/openscience/cdk/tools</i> in the folder <i>src/main/java</i> contains the class 
<i>ErtlFunctionalGroupsFinderUtilities</i> that offers static methods for required functionalities when using 
<i>ErtlFunctionalGroupsFinder</i>. It is basically a public static re-implementation of 
<i>ErtlFunctionalGroupsFinderEvaluationTest</i> from the original 
[ErtlFunctionalGroupsFinder repository](https://www.github.com/zielesny/ErtlFunctionalGroupsFinder).

The package <i>de/unijena/cheminf/ertlfgforcoconut</i> in the folder <i>src/test/java</i> contains the JUnit test class
<i>ErtlFunctionalGroupsForCoconutTest</i> that can be used to extract all Ertl functional groups from the molecules in 
COCONUT (as a MongoDB databse) and compile the information which functional groups occurred in each molecule in a text 
(CSV) file.

Other functionalities will be provided in the future.

The project is a Gradle project. In order to use it, download or clone the repository and open it in a Gradle-supporting 
IDE (e.g. IntelliJ) as a Gradle project and execute the <i>build.gradle</i> file.

## References
**ErtlFunctionalGroupsFinder**
* [Fritsch S, Neumann S, Schaub J, Steinbeck C, Zielesny A. ErtlFunctionalGroupsFinder: automated rule-based functional group detection with the Chemistry Development Kit (CDK). J Cheminform. 2019; 11:37](https://doi.org/10.1186/s13321-019-0361-8)
* [ErtlFunctionalGroupsFinder on GitHub](https://www.github.com/zielesny/ErtlFunctionalGroupsFinder)

**COCONUT**
* [Sorokina M, Steinbeck C. Review on natural products databases: where to find data in 2020. J Cheminform. 2020); 12:20](https://doi.org/10.1186/s13321-020-00424-9)
* [COCONUT on GitHub](https://github.com/mSorok/COCONUT)

**Ertl algorithm**
* [Ertl P. An algorithm to identify functional groups in organic molecules. J Cheminform. 2017; 9:36.](https://doi.org/10.1186/s13321-017-0225-z)

**Chemistry Development Kit (CDK)**
* [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* [Steinbeck C, Han Y, Kuhn S, Horlacher O, Luttmann E, Willighagen EL. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J Chem Inform Comput Sci. 2003;43(2):493-500.](https://dx.doi.org/10.1021%2Fci025584y)
* [Steinbeck C, Hoppe C, Kuhn S, Floris M, Guha R, Willighagen EL. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr Pharm Des. 2006; 12(17):2111-2120.](https://doi.org/10.2174/138161206777585274)
* [May JW and Steinbeck C. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014; 6:3.](https://dx.doi.org/10.1186%2F1758-2946-6-3)
* [Willighagen EL, Mayfield JW, Alvarsson J, Berg A, Carlsson L, Jeliazkova N, Kuhn S, Pluska T, Rojas-Chertó M, Spjuth O, Torrance G, Evelo CT, Guha R, Steinbeck C, The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J Cheminform. 2017; 9:33.](https://doi.org/10.1186/s13321-017-0220-4)
