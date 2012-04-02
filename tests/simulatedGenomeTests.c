/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
#include "sonLib.h"
#include "referenceProblem.h"

int main(int argc, char *argv[]) {
    /*
     * Parse in the chain number
     */
    int32_t nodeNumber;
    int i = scanf("%i", &nodeNumber);
    assert(i == 1);
    assert(nodeNumber >= 0 && nodeNumber % 2 == 0);
    stList *chains = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t j=0; j<nodeNumber; j++) {
        stList_append(chains, stIntTuple_construct(2, j, j+1));
    }
    /*
     * Parse in the chain weights as a list
     * of the form end1, end2, weight
     */
    double *zMatrix = st_calloc(nodeNumber*nodeNumber, sizeof(double));
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
    stList *ordering = makeReferenceGreedily(NULL, chains, zMatrix, nodeNumber, &totalScore, 0);

    /*
     * Print out the median genome,
     * as a list.
     */
}
