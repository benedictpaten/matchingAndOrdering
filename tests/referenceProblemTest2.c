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

static refAdjList *aL = NULL;
static refAdjList *dAL = NULL;
static reference *ref;
static int64_t nodeNumber;
static int64_t intervalNumber;
static int64_t weightNumber;
static int64_t testNumber = 20;

static void teardown() {
    if (aL != NULL) {
        refAdjList_destruct(aL);
        refAdjList_destruct(dAL);
        reference_destruct(ref);
        aL = NULL;
    }
}

int64_t getRandomNode(int64_t nodeNumber) {
    assert(nodeNumber > 0);
    while (1) {
        int64_t n = st_randomInt(-nodeNumber, nodeNumber + 1);
        if (n != 0) {
            return n;
        }
    }
}

static void setup() {
    teardown();
    nodeNumber = st_randomInt(0, 100) * 2;
    intervalNumber = nodeNumber > 0 ? (nodeNumber > 2 ? st_randomInt(1, nodeNumber / 2) : 1) : 0;
    weightNumber = nodeNumber >= 1 ? st_randomInt(0, nodeNumber * nodeNumber * 4) : 0;

    aL = refAdjList_construct(nodeNumber);
    dAL = refAdjList_construct(nodeNumber);
    ref = reference_construct();
    for (int64_t i = 0; i < intervalNumber; i++) {
        reference_makeNewInterval(ref, 2 * i + 1, 2 * i + 2);
    }

    for (int64_t i = 0; i < weightNumber; i++) {
        int64_t node1 = getRandomNode(nodeNumber);
        int64_t node2 = getRandomNode(nodeNumber);
        double score = st_random();
        refAdjList_addToWeight(aL, node1, node2, score);
        if (score > 0.8 && refAdjList_getWeight(dAL, node1, node2) == 0.0) {
            refAdjList_addToWeight(dAL, node1, node2, score); //0.8 thresholds ensures that nudging algorithms always results in fewer adjacencies.
        }
    }
    st_logDebug("To test the adjacency problem we've created a problem with %" PRIi64 " nodes, %" PRIi64 " intervals and %" PRIi64 " weights\n", nodeNumber,
            intervalNumber, weightNumber);
}

static void checkIsValidReference(CuTest *testCase) {
    bool *nodes = st_calloc(nodeNumber, sizeof(bool));
    for (int64_t i = 0; i < nodeNumber; i++) {
        nodes[i] = 0;
    }
    for (int64_t i = 0; i < intervalNumber; i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        CuAssertIntEquals(testCase, 2*i+1, n);
        while (reference_getNext(ref, n) != INT64_MAX) {
            CuAssertTrue(testCase, n <= nodeNumber);
            CuAssertTrue(testCase, n >= -nodeNumber);
            CuAssertTrue(testCase, n != 0);
            CuAssertTrue(testCase, nodes[abs(n)-1] == 0);
            nodes[abs(n) - 1] = 1;
            n = reference_getNext(ref, n);
        }
        CuAssertTrue(testCase, nodes[abs(n)-1] == 0);
        nodes[abs(n) - 1] = 1;
        CuAssertIntEquals(testCase, 2*i+2, n);
    }
    for (int64_t i = 0; i < nodeNumber; i++) {
        CuAssertTrue(testCase, nodes[i] == 1);
    }
    free(nodes);
}

static void testEdge(CuTest *testCase) {
    refEdge e = refEdge_construct(1, 5.0);
    CuAssertIntEquals(testCase, refEdge_to(&e), 1);
    CuAssertDblEquals(testCase, refEdge_weight(&e), 5.0, 0.0);
}

static void testAdjList(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        CuAssertIntEquals(testCase, refAdjList_getNodeNumber(aL),nodeNumber);
        if (nodeNumber > 0) {
            for (int64_t j = 0; j < 100; j++) {
                int64_t n1 = getRandomNode(nodeNumber), n2 = getRandomNode(nodeNumber);
                float w = st_random();
                refAdjList_setWeight(aL, n1, n2, w);
                CuAssertDblEquals(testCase, refAdjList_getWeight(aL, n1, n2), w, 0.0);
            }
        }
        long double totalWeight = 0.0;
        for (int64_t n1 = -nodeNumber; n1 <= nodeNumber; n1++) {
            if (n1 != 0) {
                for (int64_t n2 = n1; n2 <= nodeNumber; n2++) {
                    if (n2 != 0) {
                        double w = st_random();
                        refAdjList_setWeight(aL, n1, n2, w);
                        totalWeight += w;
                    }
                }
            }
        }
        st_logDebug("The weights were %f %f\n", refAdjList_getMaxPossibleScore(aL), totalWeight);
        CuAssertDblEquals(testCase, refAdjList_getMaxPossibleScore(aL), totalWeight, 0.0001);
        teardown();
    }
}

static void testReference(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        CuAssertIntEquals(testCase, reference_getIntervalNumber(ref), intervalNumber);
        for (int64_t j = 0; j < intervalNumber; j++) {
            CuAssertIntEquals(testCase, reference_getFirstOfInterval(ref, j), j*2 + 1);
            CuAssertIntEquals(testCase, reference_getNext(ref, j*2 + 1), j*2 + 2);
            CuAssertTrue(testCase, reference_getNext(ref, j*2 + 2) == INT64_MAX);
            CuAssertTrue(testCase, reference_getLast(ref, j*2 + 1) == j*2 + 2);
            CuAssertTrue(testCase, reference_getPrevious(ref, j*2 + 1) == INT64_MAX);
            CuAssertTrue(testCase, reference_getPrevious(ref, j*2 + 2) == j*2 + 1);
            CuAssertTrue(testCase, reference_inGraph(ref, j*2 + 1));
            CuAssertTrue(testCase, reference_inGraph(ref, j*2 + 2));
            CuAssertTrue(testCase, reference_cmp(ref, j*2 + 1, j*2 + 2) == -1);
            CuAssertTrue(testCase, reference_cmp(ref, j*2 + 2, j*2 + 1) == 1);
            CuAssertTrue(testCase, reference_cmp(ref, j*2 + 2, j*2 + 2) == 0);
            CuAssertTrue(testCase, reference_getOrientation(ref, (j*2 + 1)));
            CuAssertTrue(testCase, !reference_getOrientation(ref, -(j*2 + 1)));
            CuAssertTrue(testCase, reference_getOrientation(ref, (j*2 + 2)));
            CuAssertTrue(testCase, !reference_getOrientation(ref, -(j*2 + 2)));
            CuAssertTrue(testCase, reference_getFirst(ref, j*2 + 1) == j*2 + 1);
            CuAssertTrue(testCase, reference_getFirst(ref, j*2 + 2) == j*2 + 1);
        }
        st_logInfo("The reference for the %" PRIi64 " th test\n", i);
        reference_log(ref);
        teardown();
    }
}

static void testReferenceRandom(CuTest *testCase) {
    for (int64_t i = 0; i < testNumber; i++) {
        setup();
        time_t startTime = time(NULL);
        for (int64_t n = 2 * intervalNumber + 1; n <= nodeNumber; n++) {
            reference_insertNode(ref, n - 1 > 2 * intervalNumber ? n - 1 : 1, n);
        }
        st_logInfo("Random it took %" PRIi64 " seconds, score: %f of possible: %f\n", time(NULL) - startTime, getReferenceScore(aL, ref),
                refAdjList_getMaxPossibleScore(aL));
        st_logInfo("The reference for the %" PRIi64 " th test\n", i);
        reference_log(ref);
        checkIsValidReference(testCase);
        teardown();
    }
}

static void testMakeReferenceGreedily(CuTest *testCase) {
    long double maxScore = 0, achievedScore = 0;
    for (int64_t i = 0; i < 100; i++) {
        setup();
        time_t startTime = time(NULL);
        makeReferenceGreedily2(aL, ref, 0.99);
        int64_t badAdjacencyCount = getBadAdjacencyCount(dAL, ref);
        st_logInfo("Greedy it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n", time(NULL) - startTime,
                getReferenceScore(aL, ref), refAdjList_getMaxPossibleScore(aL), badAdjacencyCount);
        checkIsValidReference(testCase);
        updateReferenceGreedily(aL, ref, 10);
        long double greedyPermutationScore = getReferenceScore(aL, ref);
        int64_t badAdjacencyCountGreedyPermutations = getBadAdjacencyCount(dAL, ref);
        st_logInfo("Greedy with update permutations, it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n",
                time(NULL) - startTime, greedyPermutationScore, refAdjList_getMaxPossibleScore(aL), badAdjacencyCountGreedyPermutations);
        checkIsValidReference(testCase);
        reorderReferenceToAvoidBreakpoints(aL, ref);
        long double topologicalReorderedScore = getReferenceScore(aL, ref);
        checkIsValidReference(testCase);
        int64_t topologicalBadAdjacencyCount = getBadAdjacencyCount(dAL, ref);
        st_logInfo("Reordered score, it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n", time(NULL) - startTime,
                topologicalReorderedScore, refAdjList_getMaxPossibleScore(aL), topologicalBadAdjacencyCount);
        CuAssertTrue(testCase, topologicalReorderedScore >= greedyPermutationScore);
        //CuAssertTrue(testCase, getBadAdjacencyCount(dAL, ref) <= badAdjacencyCountGreedyPermutations);
        nudgeGreedily(dAL, aL, ref, 10, 100);
        long double nudgeScore = getReferenceScore(aL, ref);
        checkIsValidReference(testCase);
        st_logInfo("Nudge score, it took %" PRIi64 " seconds, score: %Lf of possible: %Lf, bad adjacency count: %" PRIi64 "\n", time(NULL) - startTime,
                nudgeScore, refAdjList_getMaxPossibleScore(aL), getBadAdjacencyCount(dAL, ref));
        CuAssertTrue(testCase, nudgeScore >= topologicalReorderedScore);
        CuAssertTrue(testCase, getBadAdjacencyCount(aL, ref) <= topologicalBadAdjacencyCount);
        reference_log(ref);
        maxScore += refAdjList_getMaxPossibleScore(aL);
        achievedScore += nudgeScore;
        teardown();
    }
    st_logInfo("Got %Lf of possible %Lf score\n", achievedScore, maxScore);
}

static void fn(double theta, int64_t node1, int64_t node2, int64_t adjacencyLength, int64_t node1Length, int64_t node2Length,
        int64_t degree) {
    double d = degree * calculateZScore(node1Length, node2Length, adjacencyLength, theta);
    assert(node1 != node2);
    refAdjList_setWeight(aL, node1, node2, refAdjList_getWeight(aL, node1, node2) + d);
}

static void testADBDCExample(CuTest *testCase) {
    /*
     * Tests example from paper.
     */
    //Nodes
    int64_t A = 1;
    int64_t AL = 2;
    int64_t C = 2;
    int64_t CL = 2;
    int64_t _5B = 3;
    int64_t _3B = -3;
    int64_t BL = 2;
    int64_t _5D = 4;
    int64_t _3D = -4;
    int64_t DL = 2;
    nodeNumber = 4;
    aL = refAdjList_construct(nodeNumber);
    ref = reference_construct();

    int64_t adjacencyLength = 1;
    int64_t n = 100;
    double theta = 0.0;

    reference_makeNewInterval(ref, C, A);
    //stList_append(chains, stIntTuple_construct(2, _5B, _3B));
    //stList_append(chains, stIntTuple_construct(2, _5D, _3D));

    fn(theta, A, _5B, 2 * adjacencyLength + DL, AL, BL, n - 1);
    fn(theta, A, _3B, 2 * adjacencyLength + DL, AL, BL, 1);
    fn(theta, A, _5D, adjacencyLength, AL, DL, n);
    fn(theta, A, _3D, 3 * adjacencyLength + DL + BL, AL, DL, n);

    fn(theta, A, C, 4 * adjacencyLength + 2 * DL + BL, AL, CL, n);

    fn(theta, C, _5B, 2 * adjacencyLength + DL, CL, BL, 1);
    fn(theta, C, _3B, 2 * adjacencyLength + DL, CL, BL, n - 1);
    fn(theta, C, _5D, adjacencyLength + DL, CL, DL, n);
    fn(theta, C, _3D, 3 * adjacencyLength + DL + BL, CL, DL, n);

    fn(theta, _3D, _5B, adjacencyLength, DL, BL, n);
    fn(theta, _3D, _3B, adjacencyLength, DL, BL, n);

    makeReferenceGreedily2(aL, ref, 0.99);
    updateReferenceGreedily(aL, ref, 100);
    st_logInfo("Running reference example problem, score: %f of possible: %f\n", getReferenceScore(aL, ref),
            refAdjList_getMaxPossibleScore(aL));
    reference_log(ref);
}

CuSuite* referenceProblem2TestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testEdge);
    SUITE_ADD_TEST(suite, testAdjList);
    SUITE_ADD_TEST(suite, testReference);
    SUITE_ADD_TEST(suite, testReferenceRandom);
    SUITE_ADD_TEST(suite, testMakeReferenceGreedily);
    SUITE_ADD_TEST(suite, testADBDCExample);
    return suite;
}
