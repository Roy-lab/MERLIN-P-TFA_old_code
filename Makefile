LFLAG = -lgsl -lgslcblas -lpthread -pthread
CC=g++
#CC=g++-4.7
CFLAGS = -g -std=c++0x
#CFLAGS = -g -std=gnu++11

OBJ = main.o Framework.o Matrix.o VariableSelection.o \
	  VariableManager.o Variable.o EvidenceManager.o Evidence.o Graph.o Error.o \
	  NCALearner.o NCA.o GraphLearner.o Potential.o PotentialManager.o SlimFactor.o \
	  FactorGraph.o FactorManager.o Distance.o HierarchicalClusterNode.o HierarchicalCluster.o \
	  LatticeStructure.o MetaMove.o HyperGeomPval.o ThreadManager.o EdgeList.o MemoryCheck.o

BIN=MERLINNCA

$(BIN): $(OBJ)
	$(CC) $(OBJ) $(LFLAG) $(CFLAGS) -o $(BIN)

main.o: main.C Framework.H
	$(CC) -c main.C $(CFLAGS) 

Framework.o: Framework.C NCALearner.H NCA.H EvidenceManager.H VariableManager.H Variable.H Graph.H EdgeList.H ThreadManager.H MemoryCheck.H
	$(CC) -c Framework.C $(CFLAGS) 

NCALearner.o: NCALearner.C NCALearner.H NCA.H EvidenceManager.H VariableManager.H Variable.H Graph.H Matrix.H ThreadManager.H
	$(CC) -c NCALearner.C $(CFLAGS) 

NCA.o: NCA.C NCA.H VariableSelection.H EvidenceManager.H VariableManager.H Variable.H Graph.H Matrix.H 
	$(CC) -c NCA.C $(CFLAGS) 

VariableSelection.o: VariableSelection.C VariableSelection.H Matrix.H
	$(CC) -c VariableSelection.C  $(CFLAGS) 

Variable.o: Variable.C Variable.H
	$(CC) -c Variable.C $(CFLAGS) 

Evidence.o: Evidence.C Evidence.H
	$(CC) -c Evidence.C $(CFLAGS) 

VariableManager.o: VariableManager.C VariableManager.H Variable.H Error.H
	$(CC) -c VariableManager.C $(CFLAGS) 

EvidenceManager.o: EvidenceManager.C EvidenceManager.H VariableManager.H Variable.H Evidence.H Matrix.H Error.H CommonTypes.H
	$(CC) -c EvidenceManager.C $(CFLAGS) 

Graph.o: Graph.C Graph.H EdgeList.H VariableManager.H Variable.H Matrix.H
	$(CC) -c Graph.C $(CFLAGS) 

Matrix.o: Matrix.C Matrix.H
	$(CC) -c Matrix.C $(CFLAGS) 

Error.o: Error.C Error.H
	$(CC) -c Error.C $(CFLAGS)

Potential.o: Potential.C Potential.H Variable.H Evidence.H
	$(CC) -c Potential.C $(CFLAGS)

PotentialManager.o: PotentialManager.C PotentialManager.H Potential.H SlimFactor.H EvidenceManager.H Variable.H  
	$(CC) -c PotentialManager.C $(CFLAGS)

SlimFactor.o: SlimFactor.C SlimFactor.H Variable.H Potential.H 
	$(CC) -c SlimFactor.C $(CFLAGS)

FactorGraph.o: FactorGraph.C FactorGraph.H EdgeList.H SlimFactor.H VariableManager.H Variable.H FactorManager.H
	$(CC) -c FactorGraph.C $(CFLAGS)

FactorManager.o: FactorManager.C FactorManager.H FactorGraph.H SlimFactor.H Error.H VariableManager.H Variable.H \
	EvidenceManager.H PotentialManager.H 
	$(CC) -c FactorManager.C $(CFLAGS)

Distance.o: Distance.C Distance.H
	$(CC) -c Distance.C $(CFLAGS)

HyperGeomPval.o: HyperGeomPval.C HyperGeomPval.H
	$(CC) -c HyperGeomPval.C $(CFLAGS)

HierarchicalClusterNode.o: HierarchicalClusterNode.C HierarchicalClusterNode.H
	$(CC) -c HierarchicalClusterNode.C $(CFLAGS)

HierarchicalCluster.o: HierarchicalCluster.C HierarchicalCluster.H HierarchicalClusterNode.H \
	VariableManager.H Variable.H Distance.H
	$(CC) -c HierarchicalCluster.C $(CFLAGS)

LatticeStructure.o: LatticeStructure.C LatticeStructure.H
	$(CC) -c LatticeStructure.C $(CFLAGS)

MetaMove.o: MetaMove.C MetaMove.H Potential.H
	$(CC) -c MetaMove.C $(CFLAGS)

GraphLearner.o: GraphLearner.C GraphLearner.H Variable.H VariableManager.H EvidenceManager.H Evidence.H \
	FactorManager.H FactorGraph.H SlimFactor.H HierarchicalCluster.H HierarchicalClusterNode.H EdgeList.H
	$(CC) -c GraphLearner.C $(CFLAGS)

EdgeList.o: EdgeList.C EdgeList.H
	$(CC) -c EdgeList.C $(CFLAGS) 

ThreadManager.o: ThreadManager.C ThreadManager.H
	$(CC) -c ThreadManager.C $(CFLAGS) 

MemoryCheck.o: MemoryCheck.C MemoryCheck.H
	$(CC) -c MemoryCheck.C $(CFLAGS) 

clean:
	rm -f $(BIN) *.o *~
