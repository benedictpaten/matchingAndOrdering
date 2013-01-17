/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stReferenceProblem2.h"
#include "stReferenceProblem.h"

static adjList *aL = NULL;
static reference *ref = NULL;
static int32_t nodeNumber = -1;
static int32_t intervalNumber = -1;
static int32_t weightNumber = -1;

static void teardown() {
    if (aL == NULL) {
        adjList_destruct(aL);
        reference_destruct(ref);
    }
}

int32_t getRandomNode(int32_t nodeNumber) {
    assert(nodeNumber > 0);
    while(1) {
        int32_t n = st_randomInt(-nodeNumber, nodeNumber + 1);
        if(nodeNumber != 0) {
            return n;
        }
    }
}

static void setup() {
    teardown();
    nodeNumber = st_randomInt(0, 100) * 2;
    intervalNumber = nodeNumber > 0 ? st_randomInt(0, nodeNumber/2) : 0;
    weightNumber = nodeNumber >= 1 ? st_randomInt(0, nodeNumber * nodeNumber * 4) : 0;

    aL = adjList_construct(nodeNumber);
    ref = reference_construct();
    for(int32_t i=0; i<intervalNumber; i++) {
        reference_makeNewInterval(ref, 2*i+1, 2*i+2);
    }

    for(int32_t i=0; i<weightNumber; i++) {
        adjList_setWeight(aL, getRandomNode(nodeNumber), getRandomNode(nodeNumber), st_random());
    }
    st_logDebug("To test the adjacency problem we've created a problem with %i nodes %i, %i intervals and %i weights\n", nodeNumber, intervalNumber, weightNumber);
}

static void checkIsValidReference(CuTest *testCase) {
    bool *nodes = st_calloc(nodeNumber, sizeof(bool));
    for(int32_t i=0; i<nodeNumber; i++) {
        nodes[i] = 0;
    }
    for(int32_t i=0; i<intervalNumber; i++) {
        int32_t n = reference_getFirstOfInterval(ref, i);
        CuAssertIntEquals(testCase, 2*i+1, n);
        while(reference_getNext(ref, n) != INT32_MAX) {
            CuAssertTrue(testCase, n <= nodeNumber);
            CuAssertTrue(testCase, n >= -nodeNumber);
            CuAssertTrue(testCase, n != 0);
            CuAssertTrue(testCase, nodes[abs(nodeNumber)-1] == 0);
            nodes[abs(nodeNumber)-1] = 1;
            n = reference_getNext(ref, n);
        }
        CuAssertTrue(testCase, nodes[abs(nodeNumber)-1] == 0);
        nodes[abs(nodeNumber)-1] = 1;
        CuAssertIntEquals(testCase, 2*i+2, n);
    }
    for(int32_t i=0; i<nodeNumber; i++) {
        CuAssertTrue(testCase, nodes[i] == 1);
    }
    free(nodes);
}

static void testMakeReferenceGreedily(CuTest *testCase) {
    for (int32_t i = 0; i < 100; i++) {
        setup();
        time_t startTime = time(NULL);
        makeReferenceGreedily2(aL, ref);
        st_logInfo("Greedy it took %i seconds, score: %f of possible: %f\n", time(NULL) - startTime, getReferenceScore(aL, ref), adjList_getMaxPossibleScore(aL));
        checkIsValidReference(testCase);
        updateReferenceGreedily(aL, ref, 100);
        st_logInfo("Greedy with update permutations, it took %i seconds, score: %f of possible: %f\n", time(NULL) - startTime, getReferenceScore(aL, ref), adjList_getMaxPossibleScore(aL));
        checkIsValidReference(testCase);
        reference_log(ref);
        teardown();
    }
}

static void fn(double theta, int32_t node1, int32_t node2,
        int32_t adjacencyLength, int32_t node1Length, int32_t node2Length, int32_t degree) {
    double d = degree * calculateZScore(node1Length, node2Length, adjacencyLength, theta);
    assert(node1 != node2);
    adjList_setWeight(aL, node1, node2, adjList_getWeight(aL, node1, node2) + d);
}

static void testADBDCExample(
        CuTest *testCase) {
    /*
     * Tests example from paper.
     */
    //Nodes
    int32_t A=1;
    int32_t AL = 2;
    int32_t C=2;
    int32_t CL = 2;
    int32_t _5B = 3;
    int32_t _3B = -3;
    int32_t BL = 2;
    int32_t _5D = 4;
    int32_t _3D = -4;
    int32_t DL = 2;
    nodeNumber = 4;
    aL = adjList_construct(nodeNumber);
    ref = reference_construct();

    int32_t adjacencyLength = 1;
    int32_t n = 100;
    double theta = 0.0;

    reference_makeNewInterval(ref, C, A);

    //stList_append(chains, stIntTuple_construct(2, _5B, _3B));
    //stList_append(chains, stIntTuple_construct(2, _5D, _3D));

    fn(theta, A, _5B, 2*adjacencyLength + DL, AL, BL, n-1);
    fn(theta, A, _3B, 2*adjacencyLength + DL, AL, BL, 1);
    fn(theta, A, _5D, adjacencyLength, AL, DL, n);
    fn(theta, A, _3D, 3*adjacencyLength + DL + BL, AL, DL, n);

    fn(theta, A, C, 4*adjacencyLength + 2*DL + BL, AL, CL, n);

    fn(theta, C, _5B, 2*adjacencyLength + DL, CL, BL, 1);
    fn(theta, C, _3B, 2*adjacencyLength + DL, CL, BL, n-1);
    fn(theta, C, _5D, adjacencyLength + DL, CL, DL, n);
    fn(theta, C, _3D, 3*adjacencyLength + DL + BL, CL, DL, n);

    fn(theta, _3D, _5B, adjacencyLength, DL, BL, n);
    fn(theta, _3D, _3B, adjacencyLength, DL, BL, n);

    makeReferenceGreedily2(aL, ref);
    updateReferenceGreedily(aL, ref, 100);
    st_logInfo("Running reference example problem, score: %f of possible: %f\n", getReferenceScore(aL, ref), adjList_getMaxPossibleScore(aL));
    reference_log(ref);
}

CuSuite* referenceProblem2TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testMakeReferenceGreedily);
    SUITE_ADD_TEST(suite, testADBDCExample);
    return suite;
}
