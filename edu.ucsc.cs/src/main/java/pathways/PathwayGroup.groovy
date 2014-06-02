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
m.add predicate:"drug" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID] 
m.add predicate:"neighborhood" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"coexpression" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"cooccurence" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"textmining" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"experimental" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"fusion" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
//m.add predicate:"database" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/* use this predicate to constrain the triad rules to only ground out valid pathways that are in pathwayPairs files*/
m.add predicate:"valPath" , types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

/*
 * Basic similarity rules for each type of similarity score
 */

initialWeight = 5

m.add rule: (neighborhood(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (coexpression(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (cooccurence(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (textmining(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (experimental(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (fusion(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
//m.add rule: (database(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (drug(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight

m.add rule: (~(neighborhood(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (~(coexpression(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (~(cooccurence(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (~(textmining(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (~(experimental(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (~(fusion(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight
//m.add rule: (database(G1, G2) & (G1 - G2)) >> pathwayPairs(G1, G2) , weight : initialWeight
m.add rule: (~(drug(G1, G2)) & (G1 - G2)) >> ~pathwayPairs(G1, G2) , weight : initialWeight



/*Triad rules for collective inference */

/*if gene a and b are linked, and c is similar to b, then a and c will be linked */
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & neighborhood(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & coexpression(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & cooccurence(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & textmining(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & experimental(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & fusion(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & database(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & drug(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight

/*if gene a and b are not linked, then a will not be linked to any c such that c is similar to b*/
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & neighborhood(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & coexpression(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & cooccurence(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & textmining(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & experimental(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & fusion(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & database(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
m.add rule: (valPath(G1, G3) & valPath(G1, G2) & ~pathwayPairs(G1, G2) & drug(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight


///*if gene a and b are linked, and c is dissimilar to b, then a will be not be linked to c */
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~neighborhood(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~coexpression(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~cooccurence(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~textmining(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~experimental(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~fusion(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
////m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & database(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~drug(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> ~pathwayPairs(G1, G3) , weight : initialWeight

/*if gene a and b are linked, and c is similar to b, then a and c will be linked */
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~neighborhood(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~coexpression(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~cooccurence(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~textmining(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~experimental(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~fusion(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
////m.add rule: (valPath(G1, G3) & pathwayPairs(G1, G2) & database(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~drug(G2, G3) & (G1 - G3) & (G2 - G3) & (G1 - G2)) >> pathwayPairs(G1, G3) , weight : initialWeight

/*
 * Transitivity
 */
 
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & pathwayPairs(G1, G2) & pathwayPairs(G2, G3) & (G1 -  G2) & (G2 - G3) & (G1 - G3)) >> pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & pathwayPairs(G2, G3) & (G1 -  G2) & (G2 - G3) & (G1 - G3)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & pathwayPairs(G1, G2) & ~pathwayPairs(G2, G3) & (G1 -  G2) & (G2 - G3) & (G1 - G3)) >> ~pathwayPairs(G1, G3) , weight : initialWeight
//m.add rule: (valPath(G1, G3) & valPath(G1, G2) & valPath(G2, G3) & ~pathwayPairs(G1, G2) & ~pathwayPairs(G2, G3) & (G1 -  G2) & (G2 - G3) & (G1 - G3)) >> pathwayPairs(G1, G3) , weight : initialWeight
/*
 * Prior that two genes are probably not in the same pathway
 */
m.add rule: (valPath(G1, G2)) >> ~(pathwayPairs(G1, G2)) , weight : initialWeight

/* Symmetry constraint on pathway groupings. I.e. pathwayPairs(G1, G2) = pathwayPairs(G2, G1) */

//m.add PredicateConstraint.Symmetric, on : pathwayPairs

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
Partition rv_te2 = new Partition(8)

//def dir = 'data'+java.io.File.separator + 'train' + java.io.File.separator;
def dir = 'data'+java.io.File.separator + 'dev' + java.io.File.separator + 'train' + java.io.File.separator;

/*
 * Load observed training data
 */
println "Loading training data..."

inserter = data.getInserter(neighborhood, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "neighborhood.txt");

inserter = data.getInserter(drug, observed_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "drug.txt");

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

//inserter = data.getInserter(database, observed_tr);
//InserterUtils.loadDelimitedDataTruth(inserter, dir + "database.txt");

inserter = data.getInserter(valPath, observed_tr);
InserterUtils.loadDelimitedData(inserter, dir + "pathwayPairs.txt");

/*
 * Load ground truth for training
 */

inserter = data.getInserter(pathwayPairs, truth_tr);
InserterUtils.loadDelimitedDataTruth(inserter, dir + "knownPathwayPairs.txt");

/*
 * Load RVs for training
 */

inserter = data.getInserter(pathwayPairs, rv_tr);
InserterUtils.loadDelimitedData(inserter, dir + "pathwayPairs.txt");

def testdir = 'data'+java.io.File.separator + 'dev' + java.io.File.separator + 'test' + java.io.File.separator;

/*
 * Load observed test data
 */
println "Loading test data..."

inserter = data.getInserter(neighborhood, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "neighborhood.txt");

inserter = data.getInserter(drug, observed_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "drug.txt");

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

//inserter = data.getInserter(database, observed_te);
//InserterUtils.loadDelimitedDataTruth(inserter, testdir + "database.txt");


inserter = data.getInserter(valPath, observed_te);
InserterUtils.loadDelimitedData(inserter, testdir + "pathwayPairs.txt");

/*
 * Load ground truth for testing
 */

inserter = data.getInserter(pathwayPairs, truth_te);
InserterUtils.loadDelimitedDataTruth(inserter, testdir + "knownPathwayPairs.txt");

/*
 * Load random variables for testing
 */

inserter = data.getInserter(pathwayPairs, rv_te);
InserterUtils.loadDelimitedData(inserter, testdir + "pathwayPairs.txt");

//inserter = data.getInserter(pathwayPairs, rv_te2);
//InserterUtils.loadDelimitedData(inserter, testdir + "unknown_pathwayPairs.txt");


/*Set up for weight learning */

Database distributionDB = data.getDatabase(predict_tr, [valPath, drug, neighborhood, coexpression, cooccurence, textmining, experimental, fusion] as Set, observed_tr);
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
Database testDB = data.getDatabase(predict_te, [valPath, drug, neighborhood, coexpression, cooccurence, textmining, experimental, fusion] as Set, observed_te);
Database testTruthDB = data.getDatabase(truth_te, [pathwayPairs] as Set)
Database rvDB = data.getDatabase(rv_te, [pathwayPairs] as Set)
//Database rvDB2 = data.getDatabase(rv_te2, [pathwayPairs] as Set)

/* This time populate testDB with all possible random variables*/
populator = new DatabasePopulator(testDB);
populator.populateFromDB(rvDB, pathwayPairs);
//populator.populateFromDB(rvDB2, pathwayPairs);

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

/* Perform discrete evaluation */
comparator = new DiscretePredictionComparator(testDB)
comparator.setBaseline(testDB)
comparator.setResultFilter(new MaxValueFilter(pathwayPairs, 1))
comparator.setThreshold(0.5) // treat best value as true as long as it is nonzero

Set<GroundAtom> groundings = Queries.getAllAtoms(testTruthDB, pathwayPairs)
int totalTestExamples = groundings.size()
DiscretePredictionStatistics stats = comparator.compare(pathwayPairs, totalTestExamples)
System.out.println("Accuracy: " + stats.getAccuracy())
System.out.println("F1: " + stats.getF1(DiscretePredictionStatistics.BinaryClass.POSITIVE))
System.out.println("Precision: " + stats.getPrecision(DiscretePredictionStatistics.BinaryClass.POSITIVE))
System.out.println("Recall: " + stats.getRecall(DiscretePredictionStatistics.BinaryClass.POSITIVE))

distributionDB.close()
labelsDB.close()
testDB.close()
testTruthDB.close()
rvDB.close()