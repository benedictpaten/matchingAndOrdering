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
    stList *stubs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t j=0; j<stubNumber; j+=2) {
        stList_append(stubs, stIntTuple_construct2(j, j+1));
    }
    stList *chains = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int64_t j=stubNumber; j<nodeNumber; j+=2) {
        stList_append(chains, stIntTuple_construct2(j, j+1));
    }
    /*
     * Parse in the chain weights as a list
     * of the form end1, end2, weight
     */
    float *zMatrix = st_calloc((nodeNumber)*(nodeNumber), sizeof(float));
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
    int64_t adjacencies[nodeNumber];
    for(int64_t j=0; j<nodeNumber; j++) {
        adjacencies[j] = INT32_MAX;
    }
    for(int64_t j=0; j<stList_length(adjacencyEdges); j++) {
        stIntTuple *edge = stList_get(adjacencyEdges, j);
        int64_t node1 = stIntTuple_getPosition(edge, 0);
        int64_t node2 = stIntTuple_getPosition(edge, 1);
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
