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
dataSet = "pathways"
ConfigManager cm = ConfigManager.getManager()
ConfigBundle cb = cm.getBundle(dataSet)

/* Uses H2 as a DataStore and stores it in a temp. directory by default */
def defaultPath = System.getProperty("java.io.tmpdir")
String dbPath = cb.getString("dbPath", defaultPath + File.separator + dataSet)
DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbPath, true), cb)

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
m.add predicate:"neighborhood" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"coexpression" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"cooccurence" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"textmining" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"experimental" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"fusion" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"database" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

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
Partition rv_tr = new Partition(6)
Partition rv_te = new Partition(7)

def dir = 'data'+java.io.File.separator + 'train' + java.io.File.separator;

/*
 * Load observed training data
 */
println "Loading training data..."

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

inserter = data.getInserter(pathwayPairs, truth_tr);
InserterUtils.loadDelimitedData(inserter, dir + "knownPathwayPairs.txt");

/*
 * Load RVs for training
 */

inserter = data.getInserter(pathwayPairs, rv_tr);
InserterUtils.loadDelimitedData(inserter, dir + "pathwayPairs.txt");

def testdir = 'data'+java.io.File.separator + 'test' + java.io.File.separator;

/*
 * Load observed test data
 */
println "Loading test data..."

inserter = data.getInserter(neighborhood, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "neighborhood.txt");

inserter = data.getInserter(coexpression, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "coexpression.txt");

inserter = data.getInserter(cooccurence, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "cooccurence.txt");

inserter = data.getInserter(textmining, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "textmining.txt");

inserter = data.getInserter(experimental, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "experimental.txt");

inserter = data.getInserter(fusion, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "fusion.txt");

inserter = data.getInserter(database, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "database.txt");

/*
 * Load ground truth for testing
 */

inserter = data.getInserter(pathwayPairs, truth_te);
InserterUtils.loadDelimitedData(inserter, testdir + "knownPathwayPairs.txt");

/*
 * Load random variables for testing
 */

inserter = data.getInserter(pathwayPairs, rv_te);
InserterUtils.loadDelimitedData(inserter, testdir + "pathwayPairs.txt");

/*Set up for weight learning */

Database distributionDB = data.getDatabase(predict_tr, [neighborhood, coexpression, cooccurence, textmining, experimental, fusion, database] as Set, observed_tr);
Database labelsDB = data.getDatabase(truth_tr, [pathwayPairs] as Set)
Database trainingRVDB = data.getDatabase(rv_tr, [pathwayPairs] as Set)

DatabasePopulator populator = new DatabasePopulator(distributionDB);
populator.populateFromDB(trainingRVDB, pathwayPairs);

println "Learning Model weights ..."

MaxLikelihoodMPE mle = new MaxLikelihoodMPE(m, distributionDB, labelsDB, cb)
mle.learn();
mle.close();

println m;

/*Set up inference*/
Database testDB = data.getDatabase(predict_te, [neighborhood, coexpression, cooccurence, textmining, experimental, fusion, database] as Set, observed_te);
Database testTruthDB = data.getDatabase(truth_te, [pathwayPairs] as Set)
Database rvDB = data.getDatabase(rv_te, [pathwayPairs] as Set)

/* This time populate testDB with all possible random variables*/
populator = new DatabasePopulator(testDB);
populator.populateFromDB(rvDB, pathwayPairs);

/* Run inference */
MPEInference mpe = new MPEInference(m, testDB, cb)
FullInferenceResult result = mpe.mpeInference()
System.out.println("Objective: " + result.getTotalWeightedIncompatibility())

/*Evaluation - compute ranking loss */
def comparator = new SimpleRankingComparator(testDB)
comparator.setBaseline(testTruthDB)

// Choosing what metrics to report
def metrics = [RankingScore.AUPRC, RankingScore.NegAUPRC,  RankingScore.AreaROC]
double [] score = new double[metrics.size()]

try {
    for (int i = 0; i < metrics.size(); i++) {
            comparator.setRankingScore(metrics.get(i))
            score[i] = comparator.compare(pathwayPairs)
    }
    //Storing the performance values of the current fold

    System.out.println("\nArea under positive-class PR curve: " + score[0])
    System.out.println("Area under negative-class PR curve: " + score[1])
    System.out.println("Area under ROC curve: " + score[2])
}
catch (ArrayIndexOutOfBoundsException e) {
    System.out.println("No evaluation data! Terminating!");
}

distributionDB.close()
labelsDB.close()
testDB.close()
testTruthDB.close()
rvDB.close()