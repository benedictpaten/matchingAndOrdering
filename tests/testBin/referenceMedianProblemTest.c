/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "sonLib.h"
#include "stReferenceProblem.h"

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
    stList *stubs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t j=0; j<stubNumber; j+=2) {
        stList_append(stubs, stIntTuple_construct(2, j, j+1));
    }
    stList *chains = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t j=stubNumber; j<nodeNumber; j+=2) {
        stList_append(chains, stIntTuple_construct(2, j, j+1));
    }
    /*
     * Parse in the chain weights as a list
     * of the form end1, end2, weight
     */
    double *zMatrix = st_calloc((nodeNumber)*(nodeNumber), sizeof(double));
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
        zMatrix[nodeNumber * node1 + node2] = weight;
        zMatrix[nodeNumber * node2 + node1] = weight;
    }
    /*
     * Compute the ordering
     */
    double totalScore;
    stList *reference = makeReferenceGreedily(stubs, chains, zMatrix, nodeNumber, &totalScore, 0);
    greedyImprovement(reference, chains, zMatrix, permutationNumber);

    /*
     * Print out the median genome,
     * as a list.
     */
    stList *adjacencyEdges = convertReferenceToAdjacencyEdges(reference);
    assert(stList_length(adjacencyEdges) == nodeNumber/2);
    int32_t adjacencies[nodeNumber];
    for(int32_t j=0; j<nodeNumber; j++) {
        adjacencies[j] = INT32_MAX;
    }
    for(int32_t j=0; j<stList_length(adjacencyEdges); j++) {
        stIntTuple *edge = stList_get(adjacencyEdges, j);
        int32_t node1 = stIntTuple_getPosition(edge, 0);
        int32_t node2 = stIntTuple_getPosition(edge, 1);
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
