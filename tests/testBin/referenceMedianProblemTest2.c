/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "sonLib.h"
#include "stReferenceProblem.h"
#include "stReferenceProblem2.h"
#include "stMatchingAlgorithms.h"

/*
 * These two functions convert from one version of the coordinates to another.
 */

int64_t convertN(int64_t n, int64_t stubNumber, int64_t nodeNumber) {
    assert(n >= 0 && n < nodeNumber);
    if(n <= stubNumber) {
        return n + 1;
    }
    if(n % 2 != 0) {
        return -(n/2 + stubNumber/2 + 1);
    }
    return n/2 + stubNumber/2 + 1;
}

int64_t convertIN(int64_t n, int64_t stubNumber, int64_t nodeNumber) {
    assert(n != 0 && n <= nodeNumber && n >= -nodeNumber);
    if(n > 0 && n <= stubNumber) {
        return n - 1;
    }
    if(n > 0) {
        return 2*n - stubNumber - 2;
    }
    return -2*n - stubNumber - 1;
}

static stList *convertReferenceToAdjacencyEdges2(reference *ref) {
    stList *edges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for(int64_t i=0; i<reference_getIntervalNumber(ref); i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        while(reference_getNext(ref, n) != INT32_MAX) {
            stList_append(edges, constructEdge(n, reference_getNext(ref, n)));
            n = -reference_getNext(ref, n);
        }
    }
    return edges;
}

int main(int argc, char *argv[]) {
    /*
     * Parse in the chain number
     */
    int64_t permutationNumber;
    int i = scanf("%" PRIi64 "", &permutationNumber);
    assert(i == 1);
    int64_t nodeNumber, stubNumber;
    i = scanf("%" PRIi64 "", &nodeNumber);
    assert(i == 1);
    i = scanf("%" PRIi64 "", &stubNumber);
    assert(i == 1);
    assert(stubNumber >= 0 && stubNumber % 2 == 0);
    assert(nodeNumber >= stubNumber && nodeNumber % 2 == 0);
    reference *ref = reference_construct();
    for(int64_t j=1; j<=stubNumber; j+=2) {
        reference_makeNewInterval(ref, j, j+1);
    }
    /*
     * Parse in the chain weights as a list
     * of the form end1, end2, weight
     */
    refAdjList *aL = adjList_construct(nodeNumber / 2 + stubNumber/2);
    int64_t weightNumber;
    i = scanf("%" PRIi64 "", &weightNumber);
    assert(i == 1);
    for(int64_t j=0; j<weightNumber; j++) {
        int64_t node1, node2;
        float weight;
        i = scanf("%" PRIi64 " %" PRIi64 " %f", &node1, &node2, &weight);
        assert(i == 3);
        assert(node1 != node2);
        assert(node1 >= 0 && node1 < nodeNumber);
        assert(node2 >= 0 && node2 < nodeNumber);
        adjList_addToWeight(aL, convertN(node1, stubNumber, nodeNumber), convertN(node2, stubNumber, nodeNumber), weight);
    }
    /*
     * Compute the ordering
     */

    makeReferenceGreedily2(aL, ref);
    updateReferenceGreedily(aL, ref, permutationNumber);

    /*
     * Print out the median genome,
     * as a list.
     */
    stList *adjacencyEdges = convertReferenceToAdjacencyEdges2(ref);
    assert(stList_length(adjacencyEdges) == nodeNumber/2);
    int64_t adjacencies[nodeNumber];
    for(int64_t j=0; j<nodeNumber; j++) {
        adjacencies[j] = INT32_MAX;
    }
    for(int64_t j=0; j<stList_length(adjacencyEdges); j++) {
        stIntTuple *edge = stList_get(adjacencyEdges, j);
        int64_t node1 = convertIN(stIntTuple_getPosition(edge, 0), stubNumber, nodeNumber);
        int64_t node2 = convertIN(stIntTuple_getPosition(edge, 1), stubNumber, nodeNumber);
        assert(node1 >= 0 && node1 <= nodeNumber);
        assert(node2 >= 0 && node2 <= nodeNumber);
        adjacencies[node1] = node2;
        adjacencies[node2] = node1;
    }
    for(int64_t j=0; j<nodeNumber; j++) {
        assert(adjacencies[j] != INT32_MAX);
        assert(adjacencies[j] >= 0 && adjacencies[j] < nodeNumber);
    }
    //Now print out the median
    int64_t nodesCovered = 0;
    for(int64_t j=0; j<stubNumber; j+=2) {
        int64_t node = adjacencies[j];
        assert(node >= 0 && node < nodeNumber);
        while(node != j+1) {
            nodesCovered+=2;
            printf("%" PRIi64 "\t", node);
            node += node % 2 == 0 ? 1 : -1;
            assert(node >= 0 && node < nodeNumber);
            node = adjacencies[node];
            assert(node >= 0 && node < nodeNumber);
            if(nodesCovered > nodeNumber) {
                assert(0);
            }
        }
    }
    printf("\n");
    assert(nodeNumber - stubNumber == nodesCovered);
}
