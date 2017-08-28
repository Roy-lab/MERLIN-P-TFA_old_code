SRC = Distance.C  Evidence.C FactorManager.C  GraphLearner.C   HyperGeomPval.C     Matrix.C       NCA.C         PotentialManager.C  Variable.C EdgeList.C  EvidenceManager.C  Framework.C      HierarchicalCluster.C      LatticeStructure.C  MemoryCheck.C  NCALearner.C  SlimFactor.C        VariableManager.C Error.C     FactorGraph.C      Graph.C          HierarchicalClusterNode.C  main.C   MetaMove.C     Potential.C   ThreadManager.C     VariableSelection.C

CC=g++
CFLAGS = -g -std=c++0x
LFLAG = -lgsl -lgslcblas -lpthread -pthread

MERLINNCA: $(SRC)
	$(CC) $(SRC) $(LFLAG) $(CFLAGS) -o MERLINNCA

clean:
	rm MERLINNCA *~
