package pathways;

//import edu.umd.cs.psl.application.inference.LazyMPEInference;
//import edu.umd.cs.psl.application.learning.weight.maxlikelihood.LazyMaxLikelihoodMPE;
import edu.umd.cs.psl.config.*
import edu.umd.cs.psl.database.DataStore
import edu.umd.cs.psl.database.Database;
import edu.umd.cs.psl.database.Partition;
import edu.umd.cs.psl.database.DatabasePopulator;
import edu.umd.cs.psl.database.ReadOnlyDatabase;
import edu.umd.cs.psl.database.rdbms.RDBMSDataStore
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver
import edu.umd.cs.psl.database.rdbms.driver.H2DatabaseDriver.Type
import edu.umd.cs.psl.model.predicate.Predicate;
import edu.umd.cs.psl.groovy.PSLModel;
import edu.umd.cs.psl.groovy.PredicateConstraint;
import edu.umd.cs.psl.groovy.SetComparison;
import edu.umd.cs.psl.model.argument.ArgumentType;
import edu.umd.cs.psl.model.argument.GroundTerm;
import edu.umd.cs.psl.model.atom.GroundAtom;
import edu.umd.cs.psl.model.function.ExternalFunction;
import edu.umd.cs.psl.ui.functions.textsimilarity.*
import edu.umd.cs.psl.ui.loading.InserterUtils;
import edu.umd.cs.psl.util.database.Queries;

import edu.umd.cs.psl.application.learning.weight.em.HardEM;

/*
 * The first thing we need to do is initialize a ConfigBundle and a DataStore
 */

ConfigManager cm = ConfigManager.getManager()
ConfigBundle config = cm.getBundle("basic-example")

println "Loading"

/* Uses H2 as a DataStore and stores it in a temp. directory by default */
String dbpath = "./test_db"
DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbpath, true), config)

PSLModel m = new PSLModel(this, data)

println "Building Model"
m.add predicate:"gene", types: [ArgumentType.UniqueID, ArgumentType.String]

m.add predicate:"pathwayMembership", types: [ArgumentType.UniqueID, ArgumentType.UniqueID]

m.add predicate:"geneCoexpression", types: [ArgumentType.UniqueID, ArgumentType.UniqueID]
m.add predicate:"drugProfileMatch", types: [ArgumentType.UniqueID, ArgumentType.UniqueID]


//If Gene A matches the drug profile of Gene B and Gene A is in Pathway P, then Gene B should be in Pathway P
m.add rule : ( drugProfileMatch(A,B) & pathwayMembership(A,P) & (A - B) ) >> pathwayMembership(B,P),  weight : 1

//If Gene A does not match the drug profile of Gene B and Gene A is in Pathway P, then Gene B should not be in Pathway P
m.add rule : ( ~drugProfileMatch(A,B) & pathwayMembership(A,P) & (A - B) ) >> ~pathwayMembership(B,P),  weight : 1


//Rule to prevent everything put into the same Pathway (ie creating a mega-pathway)


