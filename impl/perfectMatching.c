#include "sonLib.h"
#include "stSparseMatching.h"
#include "stCheckEdges.h"
#include "stMatchingAlgorithms.h"
#include "shared.h"

static void makeMatchingPerfect(stList *chosenEdges, stList *adjacencyEdges,
        stSortedSet *nodes) {
    /*
     * While the the number of edges is less than a perfect matching add random edges.
     */
    stSortedSet *attachedNodes = getNodeSetOfEdges(chosenEdges);
    stHash *nodesToAdjacencyEdges = getNodesToEdgesHash(adjacencyEdges);
    stIntTuple *pNode = NULL;
    stSortedSetIterator *it = stSortedSet_getIterator(nodes);
    stIntTuple *node;
    while((node = stSortedSet_getNext(it)) != NULL) {
        if (stSortedSet_search(attachedNodes, node) == NULL) {
            if (pNode == NULL) {
                pNode = node;
            } else {
                stList_append(chosenEdges,
                        getEdgeForNodes(stIntTuple_get(pNode, 0), stIntTuple_get(node, 0), nodesToAdjacencyEdges));
                pNode = NULL;
            }
        }
    }
    stSortedSet_destructIterator(it);
    assert(pNode == NULL);
    stSortedSet_destruct(attachedNodes);
    assert(stList_length(chosenEdges) * 2 == stSortedSet_size(nodes));
    stHash_destruct(nodesToAdjacencyEdges);
}

stList *getPerfectMatching(stSortedSet *nodes,
        stList *adjacencyEdges,
        stList *(*matchingAlgorithm)(stList *edges, int64_t nodeNumber)) {

    checkEdges(adjacencyEdges, nodes, 1, 1); //Checks edges are clique

    if (stSortedSet_size(nodes) == 0) { //Some of the following functions assume there are at least 2 nodes.
        return stList_construct();
    }

    stList *nonZeroWeightAdjacencyEdges = getEdgesWithGreaterThanZeroWeight(
                adjacencyEdges);
    stList *chosenEdges = getSparseMatching(nodes, nonZeroWeightAdjacencyEdges, matchingAlgorithm);
    stList_destruct(nonZeroWeightAdjacencyEdges);
    makeMatchingPerfect(chosenEdges, adjacencyEdges, nodes);

    st_logDebug(
                "Chosen a perfect matching with %" PRIi64 " edges, %" PRIi64 " cardinality and %" PRIi64 " weight\n",
                stList_length(chosenEdges), matchingCardinality(chosenEdges),
                matchingWeight(chosenEdges));

    return chosenEdges;
}
