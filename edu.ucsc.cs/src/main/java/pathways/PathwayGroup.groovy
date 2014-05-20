package pathways;

import java.util.Set;
import edu.umd.cs.psl.application.inference.MPEInference
import edu.umd.cs.psl.application.inference.LazyMPEInference;
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.LazyMaxLikelihoodMPE;
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.MaxLikelihoodMPE
import edu.umd.cs.psl.application.learning.weight.maxlikelihood.MaxPseudoLikelihood
import edu.umd.cs.psl.application.learning.weight.maxmargin.MaxMargin
import edu.umd.cs.psl.application.learning.weight.maxmargin.MaxMargin.NormScalingType
import edu.umd.cs.psl.application.learning.weight.random.FirstOrderMetropolisRandOM
import edu.umd.cs.psl.application.learning.weight.random.HardEMRandOM
import edu.umd.cs.psl.application.learning.weight.em.HardEM
import edu.umd.cs.psl.config.*
import edu.umd.cs.psl.core.*
import edu.umd.cs.psl.core.inference.*
import edu.umd.cs.psl.database.DataStore
import edu.umd.cs.psl.database.Database
import edu.umd.cs.psl.database.DatabasePopulator
import edu.umd.cs.psl.database.DatabaseQuery
import edu.umd.cs.psl.database.Partition
import edu.umd.cs.psl.database.ResultList
import edu.umd.cs.psl.database.rdbms.RDBMSDataStore
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver.Type
import edu.umd.cs.psl.evaluation.result.*
import edu.umd.cs.psl.evaluation.statistics.DiscretePredictionComparator
import edu.umd.cs.psl.evaluation.statistics.DiscretePredictionStatistics
import edu.umd.cs.psl.evaluation.statistics.filter.MaxValueFilter
import edu.umd.cs.psl.evaluation.statistics.RankingScore
import edu.umd.cs.psl.evaluation.statistics.SimpleRankingComparator
import edu.umd.cs.psl.groovy.*
import edu.umd.cs.psl.model.Model
import edu.umd.cs.psl.model.argument.ArgumentType
import edu.umd.cs.psl.model.argument.GroundTerm
import edu.umd.cs.psl.model.argument.UniqueID
import edu.umd.cs.psl.model.argument.Variable
import edu.umd.cs.psl.model.atom.GroundAtom
import edu.umd.cs.psl.model.atom.QueryAtom
import edu.umd.cs.psl.model.atom.RandomVariableAtom
import edu.umd.cs.psl.model.kernel.CompatibilityKernel
import edu.umd.cs.psl.model.parameters.PositiveWeight
import edu.umd.cs.psl.model.parameters.Weight
import edu.umd.cs.psl.ui.loading.*
import edu.umd.cs.psl.util.database.Queries


/*
 * The first thing we need to do is initialize a ConfigBundle and a DataStore
 */

ConfigManager cm = ConfigManager.getManager()
ConfigBundle config = cm.getBundle("pathways")

/* Uses H2 as a DataStore and stores it in a temp. directory by default */
String dbpath = "./test_db"
DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbpath, true), config)

PSLModel m = new PSLModel(this, data)

println "Building Model"

/*
 * Target predicate to predict if a gene pair belongs to the same 'cluster' or pathway bag
 */
m.add predicate:"pathwayPairs", types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/*
 * similarity scores given by the drug sensitivity matrices 
 */
m.add predicate:"drugProfileMatch", types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/*
 * similarity scores given by String-DB
 */
m.add predidcate:"neighborhood" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predidcate:"coexpression" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predidcate:"cooccurence" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predidcate:"textmining" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predidcate:"experimental" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predidcate:"fusion" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predidcate:"database" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/*
//If Gene A matches the drug profile of Gene B and Gene A is in Pathway P, then Gene B should be in Pathway P
m.add rule : ( drugProfileMatch(A,B) & pathwayMembership(A,P) & (A - B) ) >> pathwayMembership(B,P),  weight : 1

//If Gene A does not match the drug profile of Gene B and Gene A is in Pathway P, then Gene B should not be in Pathway P
m.add rule : ( ~drugProfileMatch(A,B) & pathwayMembership(A,P) & (A - B) ) >> ~pathwayMembership(B,P),  weight : 1
*/

/*
 * Basic similarity rules for each type of similarity score
 */

initialWeight = 1

m.add rule: (neighborhood(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (coexpression(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (cooccurence(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (textmining(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (experimental(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (fusion(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (database(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight

/*
 * Prior that two genes are probably not in the same pathway
 */
m.add rule: (~pathwayPairs(G1, G2)) , weight : initialWeight


/*
 * Training partitions
 */
Partition observed_tr = new Partition(0)
Partition predict_tr = new Partition(1)
Partition truth_tr = new Partition(2)

/*
 * Testing partitions
 */
Partition observed_te = new Partition(3)
Partition predict_te = new Partition(4)
Partition truth_te = new Partition(5)

/*
 * Partition for random variables
 */
Partition rvs
def dir = 'data'+java.io.File.separator;

/*
 * Load observed data
 */

inserter = data.getInserter(neighborhood, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "neighborhood.txt");

inserter = data.getInserter(coexpression, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "coexpression.txt");

inserter = data.getInserter(cooccurence, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "cooccurence.txt");

inserter = data.getInserter(textmining, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "textmining.txt");

inserter = data.getInserter(experimental, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "experimental.txt");

inserter = data.getInserter(fusion, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "fusion.txt");

inserter = data.getInserter(database, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "database.txt");

/*
 * Load ground truth for training
 */

/*
 * Load random variable data 
 */

inserter = data.getInserter(pathwayPairs, rvs);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "pathwayPairs.txt");

Database distributionDB = data.getDatabase(predict_tr, [neighborhood, coexpression, cooccurence, textmining, experimental, fusion, database], observed_tr);
Database labelsDB = data.getDatabase(truth_tr, [pathwayPairs])
Database rvDB = data.getDatabase(rvs, [pathwayPairs])

DatabasePopulator populator = new DatabasePopulator(distributionDB);
populator.populateFromDB(rvDB, pathwayPairs);

MaxLikelihoodMPE mpe = new MaxLikelihoodMPE(m, distributionDB, labelsDB, config)
mpe.learn();
mpe.close();

println m;






