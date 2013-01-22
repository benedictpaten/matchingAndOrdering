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

int32_t convertN(int32_t n, int32_t stubNumber, int32_t nodeNumber) {
    assert(n >= 0 && n < nodeNumber);
    if(n <= stubNumber) {
        return n + 1;
    }
    if(n % 2 != 0) {
        return -(n/2 + stubNumber/2 + 1);
    }
    return n/2 + stubNumber/2 + 1;
}

int32_t convertIN(int32_t n, int32_t stubNumber, int32_t nodeNumber) {
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
    for(int32_t i=0; i<reference_getIntervalNumber(ref); i++) {
        int32_t n = reference_getFirstOfInterval(ref, i);
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
    int32_t permutationNumber;
    int i = scanf("%i", &permutationNumber);
    assert(i == 1);
    int32_t nodeNumber, stubNumber;
    i = scanf("%i", &nodeNumber);
    assert(i == 1);
    i = scanf("%i", &stubNumber);
    assert(i == 1);
    assert(stubNumber >= 0 && stubNumber % 2 == 0);
    assert(nodeNumber >= stubNumber && nodeNumber % 2 == 0);
    reference *ref = reference_construct();
    for(int32_t j=1; j<=stubNumber; j+=2) {
        reference_makeNewInterval(ref, j, j+1);
    }
    /*
     * Parse in the chain weights as a list
     * of the form end1, end2, weight
     */
    adjList *aL = adjList_construct(nodeNumber / 2 + stubNumber/2);
    int32_t weightNumber;
    i = scanf("%i", &weightNumber);
    assert(i == 1);
    for(int32_t j=0; j<weightNumber; j++) {
        int32_t node1, node2;
        float weight;
        i = scanf("%i %i %f", &node1, &node2, &weight);
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
    int32_t adjacencies[nodeNumber];
    for(int32_t j=0; j<nodeNumber; j++) {
        adjacencies[j] = INT32_MAX;
    }
    for(int32_t j=0; j<stList_length(adjacencyEdges); j++) {
        stIntTuple *edge = stList_get(adjacencyEdges, j);
        int32_t node1 = convertIN(stIntTuple_getPosition(edge, 0), stubNumber, nodeNumber);
        int32_t node2 = convertIN(stIntTuple_getPosition(edge, 1), stubNumber, nodeNumber);
        assert(node1 >= 0 && node1 <= nodeNumber);
        assert(node2 >= 0 && node2 <= nodeNumber);
        adjacencies[node1] = node2;
        adjacencies[node2] = node1;
    }
    for(int32_t j=0; j<nodeNumber; j++) {
        assert(adjacencies[j] != INT32_MAX);
        assert(adjacencies[j] >= 0 && adjacencies[j] < nodeNumber);
    }
    //Now print out the median
    int32_t nodesCovered = 0;
    for(int32_t j=0; j<stubNumber; j+=2) {
        int32_t node = adjacencies[j];
        assert(node >= 0 && node < nodeNumber);
        while(node != j+1) {
            nodesCovered+=2;
            printf("%i\t", node);
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
