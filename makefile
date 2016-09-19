CC          	= g++
#STATIC	        = -static -static-libgcc -static-libstdc++
OPTIMIZE    	= -O $(STATIC)
DEBUG		= -O $(STATIC)
BioMaLL        	= ../common
BOOM 		= BOOM
GSL		= $(BOOM)/GSL
HMM		= .
OBJ	    	= obj
#CFLAGS	  	= $(OPTIMIZE) -w -fpermissive
LDFLAGS		= $(OPTIMIZE)
CFLAGS          = $(OPTIMIZE) -w -fpermissive -I$(GSLDIR)/include
LIBS	    	= -L$(BOOM) -lBOOM -lgsl -lgslcblas -lpthread 
#LIBS	    	= -L$(BOOM) -lBOOM -Wl,-Bstatic -lgsl -lgslcblas -lpthread -Wl,-Bdynamic

all: \
	BOOM \
	obj \
	parse \
	baum-welch \
	drop-continuous \
	summarize-hmm \
	train-state-labels \
	kmeans \
	get-likelihood \
	fastb-length \
	fastb-stats \
	random-HMM \
	model-combiner \
	extract-chain \
	extract-motifs \
	hmm-edit \
	hmm-extract-state \
	install-chain \
	sample \
	perl-files

.PHONY : BOOM
BOOM:
	cd BOOM ; mkdir -p obj ; make libBOOM.a

.PHONY : perl-files
perl-files:
	@chmod +x *.pl

obj:
	@mkdir obj

.PHONY : clean
clean:
	@rm -f obj/*.o BOOM/obj/*.o


$(OBJ)/model-combiner.o: \
		$(HMM)/model-combiner.C \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/model-combiner.o -c \
		$(HMM)/model-combiner.C

$(OBJ)/deterministic-trainer.o: \
		$(HMM)/deterministic-trainer.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/deterministic-trainer.o -c \
		$(HMM)/deterministic-trainer.C

$(OBJ)/baum-welch.o: \
		$(HMM)/baum-welch.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/baum-welch.o -c \
		$(HMM)/baum-welch.C

$(OBJ)/kmeans.o: \
		$(HMM)/kmeans.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/kmeans.o -c \
		$(HMM)/kmeans.C

$(OBJ)/train-state-labels.o: \
		$(HMM)/train-state-labels.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/train-state-labels.o -c \
		$(HMM)/train-state-labels.C

$(OBJ)/parallel-baum-welch.o: \
		$(HMM)/parallel-baum-welch.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/ParallelBaumWelch.H \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/parallel-baum-welch.o -c \
		$(HMM)/parallel-baum-welch.C

$(OBJ)/parallel-unscaled-BW.o: \
		$(HMM)/parallel-unscaled-BW.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/UnscaledForwardAlgorithm.H \
		$(HMM)/UnscaledBaumWelch.H \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/parallel-unscaled-BW.o -c \
		$(HMM)/parallel-unscaled-BW.C

$(OBJ)/random-HMM.o: \
		$(HMM)/random-HMM.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/Schema.H \
		$(HMM)/Schema.C \
		$(GSL)/MultiGauss.H \
		$(GSL)/MultiGauss.C \
		$(HMM)/Emission.[HC] \
		$(HMM)/EmissionSequence.[HC] \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/random-HMM.o -c \
		$(HMM)/random-HMM.C

$(OBJ)/log-likelihood.o: \
		$(HMM)/log-likelihood.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/log-likelihood.o -c \
		$(HMM)/log-likelihood.C

$(OBJ)/classify.o: \
		$(HMM)/classify.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/classify.o -c \
		$(HMM)/classify.C

$(OBJ)/hmm-substitute.o: \
		$(HMM)/hmm-substitute.C \
		$(BOOM)/Alphabet.H \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H
	g++ $(CFLAGS) -o $(OBJ)/hmm-substitute.o -c \
		$(HMM)/hmm-substitute.C

$(OBJ)/HMMbuilder.o: \
		$(HMM)/HMMbuilder.H \
		$(HMM)/HMMbuilder.C
	g++ $(CFLAGS) -o $(OBJ)/HMMbuilder.o -c \
		$(HMM)/HMMbuilder.C

$(OBJ)/HMMreader.o: \
		$(HMM)/HMMreader.H \
		$(HMM)/HMMreader.C
	g++ $(CFLAGS) -o $(OBJ)/HMMreader.o -c \
		$(HMM)/HMMreader.C

$(OBJ)/UnscaledHMMreader.o: \
		$(HMM)/UnscaledHMMreader.H \
		$(HMM)/UnscaledHMMreader.C
	g++ $(CFLAGS) -o $(OBJ)/UnscaledHMMreader.o -c \
		$(HMM)/UnscaledHMMreader.C

$(OBJ)/BackwardAlgorithm.o: \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/BackwardAlgorithm.o -c \
		$(HMM)/BackwardAlgorithm.C

$(OBJ)/UnscaledBackwardAlgorithm.o: \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H \
		$(HMM)/UnscaledBackwardAlgorithm.H \
		$(HMM)/UnscaledBackwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/UnscaledBackwardAlgorithm.o -c \
		$(HMM)/UnscaledBackwardAlgorithm.C

$(OBJ)/BaumWelchMT.o: \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/BaumWelchMT.C
	g++ $(CFLAGS) -o $(OBJ)/BaumWelchMT.o -c \
		$(HMM)/BaumWelchMT.C

$(OBJ)/ParallelBaumWelch.o: \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/ParallelBaumWelch.H \
		$(HMM)/ParallelBaumWelch.C
	g++ $(CFLAGS) -o $(OBJ)/ParallelBaumWelch.o -c \
		$(HMM)/ParallelBaumWelch.C

$(OBJ)/UnscaledBaumWelch.o: \
		$(HMM)/HMM.H \
		$(BOOM)/Sequence.H \
		$(HMM)/UnscaledForwardAlgorithm.H \
		$(HMM)/UnscaledBackwardAlgorithm.H \
		$(HMM)/UnscaledBaumWelch.H \
		$(HMM)/UnscaledBaumWelch.C
	g++ $(CFLAGS) -o $(OBJ)/UnscaledBaumWelch.o -c \
		$(HMM)/UnscaledBaumWelch.C

$(OBJ)/ForwardAlgorithm.o: \
		$(BOOM)/Sequence.H \
		$(HMM)/HMM.H \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/ForwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/ForwardAlgorithm.o -c \
		$(HMM)/ForwardAlgorithm.C

$(OBJ)/UnscaledForwardAlgorithm.o: \
		$(BOOM)/Sequence.H \
		$(HMM)/HMM.H \
		$(HMM)/UnscaledForwardAlgorithm.H \
		$(HMM)/UnscaledForwardAlgorithm.C
	g++ $(CFLAGS) -o $(OBJ)/UnscaledForwardAlgorithm.o -c \
		$(HMM)/UnscaledForwardAlgorithm.C

$(OBJ)/FastViterbi.o: \
		$(BOOM)/Sequence.H \
		$(HMM)/HMM.H \
		$(HMM)/HMMGraph.H \
		$(HMM)/FastViterbi.H \
		$(HMM)/FastViterbi.C
	g++ $(CFLAGS) -o $(OBJ)/FastViterbi.o -c \
		$(HMM)/FastViterbi.C

$(OBJ)/ComponentViterbi.o: \
		$(BOOM)/Sequence.H \
		$(HMM)/HMM.H \
		$(HMM)/HMMGraph.H \
		$(HMM)/ComponentViterbi.H \
		$(HMM)/ComponentViterbi.C
	g++ $(CFLAGS) -o $(OBJ)/ComponentViterbi.o -c \
		$(HMM)/ComponentViterbi.C

$(OBJ)/HMM.o: \
		$(HMM)/HMM.H \
		$(HMM)/HMM.C \
		$(BOOM)/Alphabet.H \
		$(BOOM)/Symbol.H
	g++ $(CFLAGS) -o $(OBJ)/HMM.o -c \
		$(HMM)/HMM.C

$(OBJ)/UnscaledHMM.o: \
		$(HMM)/UnscaledHMM.H \
		$(HMM)/UnscaledHMM.C \
		$(BOOM)/Alphabet.H \
		$(BOOM)/Symbol.H
	g++ $(CFLAGS) -o $(OBJ)/UnscaledHMM.o -c \
		$(HMM)/UnscaledHMM.C

baum-welch: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/baum-welch.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o baum-welch \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/baum-welch.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/hack-known-sites.o: \
		$(HMM)/hack-known-sites.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/hack-known-sites.o -c \
		$(HMM)/hack-known-sites.C

hack-known-sites: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/hack-known-sites.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o hack-known-sites \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/hack-known-sites.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/score-groups.o: \
		$(HMM)/score-groups.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/score-groups.o -c \
		$(HMM)/score-groups.C

score-groups: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/score-groups.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o score-groups \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/score-groups.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/drop-continuous.o: \
		$(HMM)/drop-continuous.C \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/drop-continuous.o -c \
		$(HMM)/drop-continuous.C

drop-continuous: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/drop-continuous.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o drop-continuous \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/drop-continuous.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/summarize-hmm.o: \
		$(HMM)/summarize-hmm.C \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/summarize-hmm.o -c \
		$(HMM)/summarize-hmm.C

summarize-hmm: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/summarize-hmm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o summarize-hmm \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/summarize-hmm.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

kmeans: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/kmeans.o \
		$(OBJ)/Kmeans.o
	g++ $(LDFLAGS) -o kmeans \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/kmeans.o \
		$(OBJ)/Kmeans.o \
		$(LIBS)

train-state-labels: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/train-state-labels.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o train-state-labels \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/train-state-labels.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/get-likelihood.o: \
		$(HMM)/get-likelihood.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/get-likelihood.o -c \
		$(HMM)/get-likelihood.C

get-likelihood: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/get-likelihood.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o get-likelihood \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/get-likelihood.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/fastb-length.o: \
		$(HMM)/fastb-length.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/fastb-length.o -c \
		$(HMM)/fastb-length.C

fastb-length: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/fastb-length.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o fastb-length \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/fastb-length.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/fastb-stats.o: \
		$(HMM)/fastb-stats.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/fastb-stats.o -c \
		$(HMM)/fastb-stats.C

fastb-stats: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/fastb-stats.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o fastb-stats \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/fastb-stats.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/fastb-add-wig-track.o: \
		$(HMM)/fastb-add-wig-track.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/fastb-add-wig-track.o -c \
		$(HMM)/fastb-add-wig-track.C

fastb-add-wig-track: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/fastb-add-wig-track.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o fastb-add-wig-track \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/fastb-add-wig-track.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

$(OBJ)/get-densities.o: \
		$(HMM)/get-densities.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/get-densities.o -c \
		$(HMM)/get-densities.C

get-densities: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/get-densities.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o get-densities \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/get-densities.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

parallel-baum-welch: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/parallel-baum-welch.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ParallelBaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o
	g++ $(LDFLAGS) -o parallel-baum-welch \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/parallel-baum-welch.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ParallelBaumWelch.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

parallel-unscaled-BW: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/parallel-unscaled-BW.o \
		$(OBJ)/UnscaledBackwardAlgorithm.o \
		$(OBJ)/UnscaledBaumWelch.o \
		$(OBJ)/UnscaledForwardAlgorithm.o \
		$(OBJ)/UnscaledHMM.o \
		$(OBJ)/UnscaledHMMreader.o
	g++ $(LDFLAGS) -o parallel-unscaled-BW \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/parallel-unscaled-BW.o \
		$(OBJ)/UnscaledHMMreader.o \
		$(OBJ)/UnscaledBackwardAlgorithm.o \
		$(OBJ)/UnscaledBaumWelch.o \
		$(OBJ)/UnscaledForwardAlgorithm.o \
		$(OBJ)/UnscaledHMM.o \
		$(LIBS)

random-HMM: \
		$(OBJ)/random-HMM.o \
		$(OBJ)/HMM.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o random-HMM \
		$(OBJ)/random-HMM.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(LIBS)

log-likelihood: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/log-likelihood.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o
	g++ $(LDFLAGS) -o log-likelihood \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/log-likelihood.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

classify: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/classify.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o
	g++ $(LDFLAGS) -o classify \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/classify.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

deterministic-trainer: \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/deterministic-trainer.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o
	g++ $(LDFLAGS) -o deterministic-trainer \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/deterministic-trainer.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

model-combiner: \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/model-combiner.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o
	g++ $(LDFLAGS) -o model-combiner \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/model-combiner.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/score-best-path.o:\
		$(HMM)/score-best-path.C
	$(CC) $(CFLAGS) -o $(OBJ)/score-best-path.o -c \
		$(HMM)/score-best-path.C
#--------------------------------------------------------
$(OBJ)/HMMGraph.o:\
		$(HMM)/HMMGraph.C \
		$(HMM)/HMMGraph.H
	$(CC) $(CFLAGS) -o $(OBJ)/HMMGraph.o -c \
		$(HMM)/HMMGraph.C
#--------------------------------------------------------
$(OBJ)/FastForward.o:\
		$(HMM)/FastForward.C \
		$(HMM)/FastForward.H
	$(CC) $(CFLAGS) -o $(OBJ)/FastForward.o -c \
		$(HMM)/FastForward.C
#---------------------------------------------------------
score-best-path: \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/FastForward.o \
		$(OBJ)/HMMreader.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/score-best-path.o
	$(CC) $(LDFLAGS) -o score-best-path \
		$(OBJ)/FastForward.o \
		$(OBJ)/score-best-path.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/HMMreader.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/viterbi-trainer.o:\
		viterbi-trainer.C
	$(CC) $(CFLAGS) -o $(OBJ)/viterbi-trainer.o -c \
		viterbi-trainer.C
#---------------------------------------------------------
viterbi-trainer: \
		$(OBJ)/viterbi-trainer.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/ComponentViterbi.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o
	$(CC) $(LDFLAGS) -o viterbi-trainer \
		$(OBJ)/viterbi-trainer.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/ComponentViterbi.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/parse.o:\
		parse.C
	$(CC) $(CFLAGS) -o $(OBJ)/parse.o -c \
		parse.C
#---------------------------------------------------------
parse: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/parse.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/IntervalDecoder.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o parse \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/parse.o \
		$(OBJ)/Schema.o \
		$(OBJ)/IntervalDecoder.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/test-regex.o:\
		test-regex.C
	$(CC) $(CFLAGS) -o $(OBJ)/test-regex.o -c \
		test-regex.C
#---------------------------------------------------------
test-regex: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/test-regex.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/IntervalDecoder.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o test-regex \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/test-regex.o \
		$(OBJ)/Schema.o \
		$(OBJ)/IntervalDecoder.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/split-maf.o:\
		split-maf.C
	$(CC) $(CFLAGS) -o $(OBJ)/split-maf.o -c \
		split-maf.C
#---------------------------------------------------------
split-maf: \
		$(OBJ)/split-maf.o
	$(CC) $(LDFLAGS) -o split-maf \
		$(OBJ)/split-maf.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/roc.o:\
		roc.C
	$(CC) $(CFLAGS) -o $(OBJ)/roc.o -c \
		roc.C
#---------------------------------------------------------
roc: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/roc.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o roc \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/roc.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/motif-enrichment.o:\
		motif-enrichment.C
	$(CC) $(CFLAGS) -o $(OBJ)/motif-enrichment.o -c \
		motif-enrichment.C
#---------------------------------------------------------
motif-enrichment: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/motif-enrichment.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o motif-enrichment \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/motif-enrichment.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/extract-gff-intervals.o:\
		extract-gff-intervals.C
	$(CC) $(CFLAGS) -o $(OBJ)/extract-gff-intervals.o -c \
		extract-gff-intervals.C
#---------------------------------------------------------
extract-gff-intervals: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/extract-gff-intervals.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o extract-gff-intervals \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/extract-gff-intervals.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/roc-chain.o:\
		roc-chain.C
	$(CC) $(CFLAGS) -o $(OBJ)/roc-chain.o -c \
		roc-chain.C
#---------------------------------------------------------
roc-chain: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/roc-chain.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o roc-chain \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/roc-chain.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/get-extrema.o:\
		get-extrema.C
	$(CC) $(CFLAGS) -o $(OBJ)/get-extrema.o -c \
		get-extrema.C
#---------------------------------------------------------
get-extrema: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/get-extrema.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o get-extrema \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/get-extrema.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/scan-with-chains.o:\
		scan-with-chains.C
	$(CC) $(CFLAGS) -o $(OBJ)/scan-with-chains.o -c \
		scan-with-chains.C
#---------------------------------------------------------
scan-with-chains: \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/scan-with-chains.o \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o
	$(CC) $(LDFLAGS) -o scan-with-chains \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/scan-with-chains.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/FastViterbi.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/sample.o:\
		sample.C \
		$(HMM)/HMM.H \
		$(HMM)/HMM.C
	$(CC) $(CFLAGS) -o $(OBJ)/sample.o -c \
		sample.C
#---------------------------------------------------------
sample: \
		$(OBJ)/sample.o \
		$(OBJ)/HMM.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o
	$(CC) $(LDFLAGS) -o sample \
		$(OBJ)/HMM.o \
		$(OBJ)/sample.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/Emission.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/ContinuousEmission.o:\
		ContinuousEmission.C\
		ContinuousEmission.H
	$(CC) $(CFLAGS) -o $(OBJ)/ContinuousEmission.o -c \
		ContinuousEmission.C
#---------------------------------------------------------
$(OBJ)/EmissionSequence.o:\
		EmissionSequence.C\
		EmissionSequence.H
	$(CC) $(CFLAGS) -o $(OBJ)/EmissionSequence.o -c \
		EmissionSequence.C
#---------------------------------------------------------
$(OBJ)/Emission.o:\
		Emission.C\
		Emission.H
	$(CC) $(CFLAGS) -o $(OBJ)/Emission.o -c \
		Emission.C
#---------------------------------------------------------

$(OBJ)/EmissionLoader.o:\
		EmissionLoader.C\
		EmissionLoader.H
	$(CC) $(CFLAGS) -o $(OBJ)/EmissionLoader.o -c \
		EmissionLoader.C
#---------------------------------------------------------
$(OBJ)/test-multigauss.o:\
		test-multigauss.C
	$(CC) $(CFLAGS) -o $(OBJ)/test-multigauss.o -c \
		test-multigauss.C
#---------------------------------------------------------
test-multigauss: \
		$(OBJ)/test-multigauss.o
	$(CC) $(LDFLAGS) -o test-multigauss \
		$(OBJ)/test-multigauss.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/Schema.o:\
		Schema.C\
		Schema.H
	$(CC) $(CFLAGS) -o $(OBJ)/Schema.o -c \
		Schema.C
#---------------------------------------------------------
$(OBJ)/Kmeans.o:\
		Kmeans.C\
		Kmeans.H
	$(CC) $(CFLAGS) -o $(OBJ)/Kmeans.o -c \
		Kmeans.C
#---------------------------------------------------------
$(OBJ)/cluster-fastb.o:\
		cluster-fastb.C
	$(CC) $(CFLAGS) -o $(OBJ)/cluster-fastb.o -c \
		cluster-fastb.C
#---------------------------------------------------------
cluster-fastb: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(OBJ)/cluster-fastb.o
	$(CC) $(LDFLAGS) -o cluster-fastb \
		$(OBJ)/cluster-fastb.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/scale-data.o:\
		scale-data.C
	$(CC) $(CFLAGS) -o $(OBJ)/scale-data.o -c \
		scale-data.C
#---------------------------------------------------------
scale-data: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/scale-data.o
	$(CC) $(LDFLAGS) -o scale-data \
		$(OBJ)/scale-data.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/Schema.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/discretize.o:\
		discretize.C
	$(CC) $(CFLAGS) -o $(OBJ)/discretize.o -c \
		discretize.C
#---------------------------------------------------------
discretize: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/discretize.o
	$(CC) $(LDFLAGS) -o discretize \
		$(OBJ)/discretize.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/Schema.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/threshold.o:\
		threshold.C
	$(CC) $(CFLAGS) -o $(OBJ)/threshold.o -c \
		threshold.C
#---------------------------------------------------------
threshold: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/Schema.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/threshold.o
	$(CC) $(LDFLAGS) -o threshold \
		$(OBJ)/threshold.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/Schema.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/SequenceSet.o:\
		SequenceSet.C\
		SequenceSet.H
	$(CC) $(CFLAGS) -o $(OBJ)/SequenceSet.o -c \
		SequenceSet.C
#---------------------------------------------------------
$(OBJ)/smooth-fastb.o:\
		smooth-fastb.C
	$(CC) $(CFLAGS) -o $(OBJ)/smooth-fastb.o -c \
		smooth-fastb.C
#---------------------------------------------------------
smooth-fastb: \
		$(OBJ)/smooth-fastb.o
	$(CC) $(LDFLAGS) -o smooth-fastb \
		$(OBJ)/smooth-fastb.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/subset-fastb.o:\
		subset-fastb.C
	$(CC) $(CFLAGS) -o $(OBJ)/subset-fastb.o -c \
		subset-fastb.C
#---------------------------------------------------------
subset-fastb: \
		$(OBJ)/subset-fastb.o
	$(CC) $(LDFLAGS) -o subset-fastb \
		$(OBJ)/subset-fastb.o \
		$(LIBS)
#---------------------------------------------
$(OBJ)/extract-motifs.o: \
		$(HMM)/extract-motifs.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/extract-motifs.o -c \
		$(HMM)/extract-motifs.C
extract-motifs: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/extract-motifs.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o extract-motifs \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/extract-motifs.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)
#---------------------------------------------
$(OBJ)/extract-chain.o: \
		$(HMM)/extract-chain.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/extract-chain.o -c \
		$(HMM)/extract-chain.C
extract-chain: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/extract-chain.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o extract-chain \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/extract-chain.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)
#---------------------------------------------
$(OBJ)/install-chain.o: \
		$(HMM)/install-chain.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/install-chain.o -c \
		$(HMM)/install-chain.C
install-chain: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/install-chain.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o install-chain \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/install-chain.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)
#---------------------------------------------
$(OBJ)/extract-state-sequences.o: \
		$(HMM)/extract-state-sequences.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/BaumWelchMT.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/extract-state-sequences.o -c \
		$(HMM)/extract-state-sequences.C
extract-state-sequences: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/extract-state-sequences.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o extract-state-sequences \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/extract-state-sequences.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)
#---------------------------------------------
$(OBJ)/entropy.o: \
		$(HMM)/entropy.C
	g++ $(CFLAGS) -o $(OBJ)/entropy.o -c \
		$(HMM)/entropy.C
entropy: \
		$(OBJ)/entropy.o
	g++ $(LDFLAGS) -o entropy \
		$(OBJ)/entropy.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/TieProfile.o:\
		TieProfile.C\
		TieProfile.H
	$(CC) $(CFLAGS) -o $(OBJ)/TieProfile.o -c \
		TieProfile.C
#---------------------------------------------------------
$(OBJ)/TieEvent.o:\
		TieEvent.C\
		TieEvent.H
	$(CC) $(CFLAGS) -o $(OBJ)/TieEvent.o -c \
		TieEvent.C
#---------------------------------------------------------

$(OBJ)/hmm-edit.o: \
		$(HMM)/hmm-edit.C \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/hmm-edit.o -c \
		$(HMM)/hmm-edit.C

hmm-edit: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/hmm-edit.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o hmm-edit \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/hmm-edit.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

#---------------------------------------------------------

$(OBJ)/hmm-extract-state.o: \
		$(HMM)/hmm-extract-state.C \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/hmm-extract-state.o -c \
		$(HMM)/hmm-extract-state.C

hmm-extract-state: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/hmm-extract-state.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o hmm-extract-state \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/SequenceSet.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/TieProfile.o \
		$(OBJ)/TieEvent.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/hmm-extract-state.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/BaumWelchMT.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)


$(OBJ)/find-peaks.o: \
		$(HMM)/find-peaks.C \
		$(HMM)/ForwardAlgorithm.H \
		$(HMM)/BackwardAlgorithm.H \
		$(HMM)/HMM.H
	g++ $(CFLAGS) -o $(OBJ)/find-peaks.o -c \
		$(HMM)/find-peaks.C

find-peaks: \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/find-peaks.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(OBJ)/HMMbuilder.o
	g++ $(LDFLAGS) -o find-peaks \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/Schema.o \
		$(OBJ)/HMMGraph.o \
		$(OBJ)/find-peaks.o \
		$(OBJ)/HMMbuilder.o \
		$(OBJ)/BackwardAlgorithm.o \
		$(OBJ)/ForwardAlgorithm.o \
		$(OBJ)/HMM.o \
		$(LIBS)

#--------------------------------------------------------
$(OBJ)/IntervalDecoder.o:\
		IntervalDecoder.C\
		IntervalDecoder.H
	$(CC) $(CFLAGS) -o $(OBJ)/IntervalDecoder.o -c \
		IntervalDecoder.C
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/fastb-slice-by-symbol.o:\
		fastb-slice-by-symbol.C
	$(CC) $(CFLAGS) -o $(OBJ)/fastb-slice-by-symbol.o -c \
		fastb-slice-by-symbol.C
#---------------------------------------------------------
fastb-slice-by-symbol: \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/fastb-slice-by-symbol.o
	$(CC) $(LDFLAGS) -o fastb-slice-by-symbol \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/fastb-slice-by-symbol.o \
		$(LIBS)
#--------------------------------------------------------
$(OBJ)/fastb-slice.o:\
		fastb-slice.C
	$(CC) $(CFLAGS) -o $(OBJ)/fastb-slice.o -c \
		fastb-slice.C
#---------------------------------------------------------
fastb-slice: \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/fastb-slice.o
	$(CC) $(LDFLAGS) -o fastb-slice-by-symbol \
		$(OBJ)/Schema.o \
		$(OBJ)/Emission.o \
		$(OBJ)/EmissionSequence.o \
		$(OBJ)/EmissionLoader.o \
		$(OBJ)/fastb-slice.o \
		$(LIBS)
#---------------------------------------------


