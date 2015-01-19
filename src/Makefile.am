# --enable-profile
# --enable-debug
# --enable-test

## i386 and x86-64 options: -minline-all-stringops -momit-leaf-frame-pointer

if DEBUG
OPTIMIZATIONS = -ggdb -DDEBUG
else
if PROFILE
OPTIMIZATIONS = -ggdb -pg -DNDEBUG -O2 -funroll-all-loops $(ARCH)
else
if TESTING
OPTIMIZATIONS = -ggdb
else
OPTIMIZATIONS = -DNDEBUG -O2 -funroll-all-loops $(ARCH)
endif
endif
endif


AM_CFLAGS = -W -Wformat -Wimplicit -Wreturn-type -Wall -Werror \
            -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
            -Winline -Wnested-externs -Wextra $(OPTIMIZATIONS) \
            -std=c99

################################################################################

bin_PROGRAMS = assemble assemble1 assemble_illumina build_sff_index \
               cull_components filter filter_clusters fish_clusters \
               make_template masker remove_clusters update update_clusters

################################################################################

CLUSTER_SRC = cluster.c cluster.h 
COMMON_SRC = common.c common.h
CONTIG_SRC = contig.c contig.h
GRAPH_SRC = graph.c graph.h
HASHTABLE_SRC = hashtable.c hashtable.h
MAF_SRC = maf.c maf.h mfblocks.c mfblocks.h
SEQUENCES_SRC = sequences.c sequences.h
SFF_SRC = sff.c sff.h
SFF_INDEX_SRC = sff.c sff.h sff_index.c sff_index.h
SLINKLIST_SRC = slinklist.c slinklist.h
SPLAYHASH_SRC = splayhash.c splayhash.h
UTILITIES_SRC = utilities.c utilities.h

LIB_NOSEQ_SRC = $(UTILITIES_SRC) $(COMMON_SRC) $(SLINKLIST_SRC) $(HASHTABLE_SRC)
LIB_SRC = $(LIB_NOSEQ_SRC) $(SEQUENCES_SRC)

################################################################################

assemble_SOURCES = $(LIB_SRC) $(SFF_SRC) $(CONTIG_SRC) assemble.c

assemble1_SOURCES = $(LIB_SRC) $(SFF_INDEX_SRC) $(CONTIG_SRC) assemble1.c

assemble_illumina_SOURCES = $(LIB_SRC) assemble_illumina.c

build_sff_index_SOURCES = $(LIB_NOSEQ_SRC) $(SFF_INDEX_SRC) build_sff_index.c

cull_components_SOURCES = $(LIB_SRC) $(GRAPH_SRC) $(CLUSTER_SRC) \
                          cull_components.c

filter_SOURCES = $(UTILITIES_SRC) filter.c

filter_clusters_SOURCES = $(LIB_NOSEQ_SRC) $(CLUSTER_SRC) filter_clusters.c

fish_clusters_SOURCES = $(LIB_SRC) $(GRAPH_SRC) $(CLUSTER_SRC) fish_clusters.c

make_template_SOURCES = $(UTILITIES_SRC) $(COMMON_SRC) $(SEQUENCES_SRC) \
                        make_template.c

masker_SOURCES = $(UTILITIES_SRC) $(COMMON_SRC) $(SPLAYHASH_SRC) \
                 $(SEQUENCES_SRC) masker.c
masker_LDADD = -lm

remove_clusters_SOURCES = $(UTILITIES_SRC) $(CLUSTER_SRC) $(SLINKLIST_SRC) \
                          $(MAF_SRC) remove_clusters.c

update_SOURCES = $(LIB_SRC) update.c

update_clusters_SOURCES = $(LIB_SRC) $(MAF_SRC) $(CLUSTER_SRC) update_clusters.c
update_clusters_LDADD = -lm

################################################################################

bin_SCRIPTS = DIAL convertFastqFasta partition remove_clones update_status

################################################################################
