/*
 * referenceProblem2.c
 *
 *  Created on: 15 Jan 2013
 *      Author: benedictpaten
 */

#include "stReferenceProblem2.h"
#include "sonLib.h"
#include <stdlib.h>
#include <math.h>

/*
 * Adjacency list structure/edge structure
 */

struct _edge {
    int32_t to;
    float weight;
};

edge edge_construct(int32_t to, float weight) {
    edge e;
    e.to = to;
    e.weight = weight;
    return e;
}

int32_t edge_to(edge *e) {
    return e->to;
}

float edge_weight(edge *e) {
    return e->weight;
}

static edge *edge_copy(edge *e) {
    edge *e2 = st_malloc(sizeof(edge));
    *e2 = *e;
    return e2;
}

static int edge_cmp(edge *e1, edge *e2, reference *ref) {
    return reference_cmp(ref, edge_to(e1), edge_to(e2));
}

struct _adjList {
    stHash **edgeHashes;
    int32_t nodeNumber;
};

adjList *adjList_construct(int32_t nodeNumber) {
    adjList *aL = st_malloc(sizeof(adjList));
    aL->nodeNumber = nodeNumber;
    aL->edgeHashes = st_malloc(sizeof(stHash *) * nodeNumber * 2);
    for (int32_t i = 0; i < 2 * nodeNumber; i++) {
        aL->edgeHashes[i] = stHash_construct3((uint32_t(*)(const void *)) stIntTuple_hashKey,
                (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct, free);
    }
    return aL;
}

void adjList_destruct(adjList *aL) {
    for (int32_t i = 0; i < 2 * aL->nodeNumber; i++) {
        stHash_destruct(aL->edgeHashes[i]);
    }
    free(aL);
}

int32_t adjList_getNodeNumber(adjList *aL) {
    return aL->nodeNumber;
}

static void checkN(int32_t n, int32_t nodeNumber) {
    assert(n <= aL->nodeNumber);
    assert(n != 0);
    assert(n >= -aL->nodeNumber);
}

static int32_t convertN(adjList *aL, int32_t n) {
    checkN(n, aL->nodeNumber);
    return (n < 0 ? -n : n + aL->nodeNumber) - 1;
}

float adjList_getWeight(adjList *aL, int32_t n1, int32_t n2) {
    checkN(n2, aL->nodeNumber);
    stIntTuple *i = stIntTuple_construct(1, n2);
    float *weight = stHash_search(aL->edgeHashes[convertN(aL, n1)], i);
    stIntTuple_destruct(i);
    return weight == NULL ? 0.0 : weight[0];
}

static void adjList_setWeightP(adjList *aL, int32_t n1, int32_t n2, float weight) {
    stHash *edges = aL->edgeHashes[convertN(aL, n1)];
    checkN(n2, aL->nodeNumber);
    stIntTuple *i = stIntTuple_construct(1, n2);
    float *w = stHash_search(edges, i);
    if (w == NULL) {
        float *w = st_malloc(sizeof(float));
        w[0] = weight;
        stHash_insert(edges, i, w);
    } else {
        w[0] = weight;
        stIntTuple_destruct(i);
    }
}

void adjList_setWeight(adjList *aL, int32_t n1, int32_t n2, float weight) {
    adjList_setWeightP(aL, n1, n2, weight);
    adjList_setWeightP(aL, n2, n1, weight);
}

struct _adjListIt {
    stHash *hash;
    stHashIterator *it;
};

adjListIt adjList_getEdgeIt(adjList *aL, int32_t node) {
    adjListIt it;
    it.hash = aL->edgeHashes[convertN(aL, node)];
    it.it = stHash_getIterator(it.hash);
    return it;
}

edge adjListIt_getNext(adjListIt *it) {
    edge e;
    stIntTuple *i = stHash_getNext(it->it);
    if (i == NULL) {
        e.to = INT32_MAX;
    } else {
        e.to = stIntTuple_getPosition(i, 0);
        e.weight = *(float *) stHash_search(it->hash, i);
    }
    return e;
}

void adjListIt_destruct(adjListIt *it) {
    stHash_destructIterator(it->it);
}

double adjList_getMaxPossibleScoreP(adjList *aL, int32_t n) {
    double score = 0.0;
    adjListIt it = adjList_getEdgeIt(aL, n);
    edge e = adjListIt_getNext(&it);
    while (edge_to(&e) != INT32_MAX) {
        score += edge_weight(&e);
    }
    adjListIt_destruct(&it);
    return score;
}

double adjList_getMaxPossibleScore(adjList *aL) {
    double score = 0.0;
    for (int32_t n = 1; n <= adjList_getNodeNumber(aL); n++) {
        score += adjList_getMaxPossibleScoreP(aL, n);
        score += adjList_getMaxPossibleScoreP(aL, -n);
    }
    return score / 2.0;
}

/*double calculateZScore(int32_t n, int32_t m, int32_t k, double theta) {
    assert(theta <= 1.0);
    assert(theta >= 0.0);
    if (theta == 0.0) {
        return ((double) n) * m;
    }
    double beta = 1.0 - theta;
    return ((1.0 - pow(beta, n)) / theta) * pow(beta, k) * ((1.0 - pow(beta, m)) / theta);
}*/

/*
 * Insertion point
 */

struct _insertPoint {
    int32_t node;
    int32_t adjNode;
    bool left;
    float score;
};

typedef struct _insertPoint insertPoint;

static insertPoint *insertPoint_construct(int32_t node, int32_t adjNode, bool left, float score) {
    insertPoint *insertPoint = st_malloc(sizeof(insertPoint));
    insertPoint->node = node;
    insertPoint->adjNode = adjNode;
    insertPoint->left = left;
    insertPoint->score = score;
    return insertPoint;
}

static int32_t insertPoint_node(insertPoint *iP) {
    return iP->node;
}

static int32_t insertPoint_adjNode(insertPoint *iP) {
    return iP->adjNode;
}

static bool insertPoint_left(insertPoint *iP) {
    return iP->left;
}

static float insertPoint_score(insertPoint *iP) {
    return iP->score;
}

static int insertPoint_cmp(insertPoint *iP1, insertPoint *iP2, reference *ref) {
    return reference_cmp(ref, insertPoint_adjNode(iP1), insertPoint_adjNode(iP2));
}

/*
 * Reference structure
 */

typedef struct _referenceTerm referenceTerm;

struct _referenceTerm {
    referenceTerm *first, *pTerm, *nTerm;
    int32_t node;
    int64_t index;
};

struct _reference {
    stHash *nodesInGraph;
    stList *referenceIntervals;
};

reference *reference_construct() {
    reference *ref = st_malloc(sizeof(reference));
    ref->nodesInGraph = stHash_construct();
    ref->referenceIntervals = stList_construct();
    return ref;
}

void reference_destruct(reference *ref) {
    stHash_destruct(ref->nodesInGraph);
    stList_destruct(ref->referenceIntervals);
    free(ref);
}

void reference_makeNewInterval(reference *ref, int32_t leftNode, int32_t rightNode) {
    referenceTerm *rTL = st_malloc(sizeof(referenceTerm)), *rTR = st_malloc(sizeof(referenceTerm));
    rTL->node = leftNode;
    rTR->node = rightNode;
    rTL->nTerm = rTR;
    rTR->pTerm = rTL;
    rTL->pTerm = NULL;
    rTR->nTerm = NULL;
    rTL->first = rTL;
    rTR->first = rTL;
    rTL->index = 0;
    rTR->index = INT64_MAX;
    stHash_insert(ref->nodesInGraph, rTL, rTL);
    stHash_insert(ref->nodesInGraph, rTR, rTR);
}

static referenceTerm *reference_getTerm(reference *ref, int32_t n) {
    referenceTerm rT;
    rT.node = n;
    return stHash_search(ref->nodesInGraph, &rT);
}

static void reference_insertNodeP(reference *ref, int32_t pNode, int32_t node) {
    referenceTerm *rTR = st_malloc(sizeof(referenceTerm)), *rTL;
    rTR->node = node;
    rTL = reference_getTerm(ref, pNode);
    assert(rTL != NULL);
    rTR->nTerm = rTL->nTerm;
    rTR->pTerm = rTL;
    rTL->nTerm = rTR;
    rTR->nTerm->pTerm = rTR;
    rTR->first = rTL->first;
    stHash_insert(ref->nodesInGraph, rTR, rTR);
    //Deal with indices
    if(rTL->nTerm->index - rTL->index == 1) { //Need to rebalance
        //Work out the length of the chain
        referenceTerm *rT = rTL->first;
        int32_t length = 0;
        while(rT != NULL) {
            length++;
            rT = rT->nTerm;
        }
        //Now give every one equi-distant labels.
        int64_t spacer = INT64_MAX/length;
        rT = rTL->first;
        rT->index = 0;
        while (rT->nTerm != NULL) {
            rT->nTerm->index = rT->index + spacer;
            rT = rT->nTerm;
        }
    }
    assert(rTL->nTerm->index - rTL->index > 1);
    rTR->index = rTL->nTerm->index + (rTL->nTerm->index - rTL->index)/2;
}

static void reference_insertNode(reference *ref, insertPoint *iP) {
    if (insertPoint_left(iP)) {
        reference_insertNodeP(ref, insertPoint_adjNode(iP), insertPoint_node(iP));
    } else {
        reference_insertNodeP(ref, reference_getPrevious(ref, insertPoint_adjNode(iP)), insertPoint_node(iP));
    }
}

static void reference_removeNode(reference *ref, int32_t n) {
    referenceTerm *rT = stHash_remove(ref->nodesInGraph, reference_getTerm(ref, n));
    rT->nTerm->pTerm = rT->pTerm;
    rT->pTerm->nTerm = rT->nTerm;
}

bool reference_inGraph(reference *ref, int32_t n) {
    return reference_getTerm(ref, n) != NULL;
}

int32_t reference_getFirstOfInterval(reference *ref, int32_t interval) {
    return ((referenceTerm *)stList_get(ref->referenceIntervals, interval))->node;
}

int32_t reference_getIntervalNumber(reference *ref) {
    return stList_length(ref->referenceIntervals);
}

int32_t reference_getFirst(reference *ref, int32_t n) {
    return reference_getTerm(ref, n)->first->node;
}

int32_t reference_getPrevious(reference *ref, int32_t n) {
    referenceTerm *rT = reference_getTerm(ref, n);
    return rT->pTerm != NULL ? rT->pTerm->node : INT32_MAX;
}

bool reference_getOrientation(reference *ref, int32_t n) {
    return reference_getTerm(ref, n)->node == n;
}

int32_t reference_getNext(reference *ref, int32_t n) {
    referenceTerm *rT = reference_getTerm(ref, n);
    return rT->nTerm != NULL ? rT->nTerm->node : INT32_MAX;
}

int reference_cmp(reference *ref, int32_t n1, int32_t n2) {
    referenceTerm *rT1 = reference_getTerm(ref, n1), *rT2 = reference_getTerm(ref, n2);
    if (rT1->first != rT2->first) {
        return rT1->first > rT2->first ? 1 : -1;
    }
    return rT1->index > rT2->index ? 1 : rT1->index < rT2->index ? -1 : 0;
}

void reference_log(reference *ref) {
    st_logInfo("Logging reference with %i intervals\n", reference_getIntervalNumber(ref));
    for(int32_t i=0; i<reference_getIntervalNumber(ref); i++) {
        st_logInfo("Interval %i, nodes:", i);
        int32_t n = reference_getFirstOfInterval(ref, i);
        while(n != INT32_MAX) {
            st_logInfo(" %i, ", n);
            n = reference_getNext(ref, n);
        }
        st_logInfo("\n");
    }
}

/*
 * Reference algorithm
 */

static stList *getRelevantEdges(adjList *aL, reference *ref, int32_t n) {
    stList *edges = stList_construct3(0, free);
    adjListIt it = adjList_getEdgeIt(aL, n);
    edge e = adjListIt_getNext(&it);
    while (edge_to(&e) != INT32_MAX) {
        //Check if edge is in graph
        if (reference_inGraph(ref, edge_to(&e))) {
            //Is in graph, so add it to list
            stList_append(edges, edge_copy(&e));
        }
        e = adjListIt_getNext(&it);
    }
    adjListIt_destruct(&it);
    //Now do sorting to determine ordering
    stList_sort2(edges, (int (*)(const void *, const void *, const void *))edge_cmp, ref);
    return edges;
}

static void getInsertPointsRight(int32_t n, stList *edges, reference *ref, stList *insertPoints) {
    double f = 0.0;
    int32_t intervalName = INT32_MAX;
    for (int32_t i = 0; i < stList_length(edges); i++) {
        edge *e = stList_get(edges, i);
        if (reference_getFirst(ref, edge_to(e)) != intervalName) { //Reset placement
            f = 0.0;
            intervalName = reference_getFirst(ref, edge_to(e));
        }
        if (!reference_getOrientation(ref, edge_to(e))) {
            f += edge_weight(e);
            stList_append(insertPoints, insertPoint_construct(n, edge_to(e), 0, f));
        }
    }
}

static void getInsertPointsLeft(int32_t n, stList *edges, reference *ref, stList *insertPoints) {
    double f = 0.0;
    int32_t intervalName = INT32_MAX;
    for (int32_t i = stList_length(edges) - 1; i >= 0; i--) {
        edge *e = stList_get(edges, i);
        if (reference_getFirst(ref, edge_to(e)) != intervalName) { //Reset placement
            f = 0.0;
            intervalName = reference_getFirst(ref, edge_to(e));
        }
        if (reference_getOrientation(ref, edge_to(e))) {
            f += edge_weight(e);
            stList_append(insertPoints, insertPoint_construct(n, edge_to(e), 1, f));
        }
    }
}

static void getInsertionPoints(int32_t n, stList *leftEdges, stList *rightEdges, reference *ref, adjList *aL, stList *insertPoints) {
    stList *insertPoints2 = stList_construct();
    getInsertPointsRight(n, leftEdges, ref, insertPoints2);
    getInsertPointsLeft(n, rightEdges, ref, insertPoints2);
    stList_sort2(insertPoints2, (int (*)(const void *, const void *, const void *))insertPoint_cmp, ref);
    int32_t j = stList_length(insertPoints2);
    for (int32_t i = 1; i < j; i++) {
        insertPoint *iPL = stList_get(insertPoints2, i - 1);
        insertPoint *iPR = stList_get(insertPoints2, i);
        if (!insertPoint_left(iPL) && insertPoint_left(iPR) && reference_getFirst(ref,
                insertPoint_adjNode(iPL) == reference_getFirst(ref, insertPoint_adjNode(iPR)))) {
            if (adjList_getWeight(aL, n, insertPoint_adjNode(iPL)) > adjList_getWeight(aL, -n, insertPoint_adjNode(iPR))) {
                stList_append(insertPoints,
                        insertPoint_construct(n, insertPoint_adjNode(iPL), 0, insertPoint_score(iPL) + insertPoint_score(iPR)));
            } else {
                stList_append(insertPoints,
                        insertPoint_construct(n, insertPoint_adjNode(iPR), 1, insertPoint_score(iPL) + insertPoint_score(iPR)));
            }
        }
    }
    stList_appendAll(insertPoints, insertPoints2);
    stList_destruct(insertPoints2);
}

static void insertNode(int32_t n, adjList *aL, reference *ref) {
    if (!reference_inGraph(ref, n)) { //Have a node to add
        //Get list of edges to nodes already in the reference
        stList *fEdges = getRelevantEdges(aL, ref, n);
        stList *nEdges = getRelevantEdges(aL, ref, -n);
        //Now do the hardwork of determining the best insertion point
        stList *insertPoints = stList_construct3(0, free);
        getInsertionPoints(n, fEdges, nEdges, ref, aL, insertPoints);
        getInsertionPoints(-n, nEdges, fEdges, ref, aL, insertPoints);
        //Get best insertion point
        insertPoint *bestIP = NULL;
        for (int32_t i = 0; i < stList_length(insertPoints); i++) {
            insertPoint *iP = stList_get(insertPoints, i);
            if (bestIP == NULL || insertPoint_score(iP) > insertPoint_score(bestIP)) {
                bestIP = iP;
            }
        }
        assert(bestIP != NULL);
        reference_insertNode(ref, bestIP);
        //Cleanup
        stList_destruct(fEdges);
        stList_destruct(nEdges);
        stList_destruct(insertPoints);
    }
}

void makeReferenceGreedily2(adjList *aL, reference *ref) {
    for (int32_t n; n <= adjList_getNodeNumber(aL); n++) {
        insertNode(n, aL, ref);
    }
}

void updateReferenceGreedily(adjList *aL, reference *ref, int32_t permutations) {
    for (int32_t i = 0; i < permutations; i++) {
        for (int32_t j; j <= adjList_getNodeNumber(aL); j++) {
            int32_t n = st_randomInt(0, adjList_getNodeNumber(aL));
            assert(reference_inGraph(ref, n));
            if (abs(reference_getFirst(ref, n)) != n) {
                reference_removeNode(ref, n);
                insertNode(n, aL, ref);
            }
        }
    }
}

static double getSumOfConsistentAdjacenciesScore(int32_t n, adjList *aL, reference *ref) {
    double score = 0.0;
    adjListIt it = adjList_getEdgeIt(aL, n);
    edge e = adjListIt_getNext(&it);
    while (edge_to(&e) != INT32_MAX) {
        if (reference_cmp(ref, n, edge_to(&e)) == -1 && !reference_getOrientation(ref, n) && reference_getOrientation(ref, edge_to(&e))) {
            score += edge_weight(&e);
        }
        e = adjListIt_getNext(&it);
    }
    adjListIt_destruct(&it);
    return score;
}

double getReferenceScore(adjList *aL, reference *ref) {
    double score = 0.0;
    for (int32_t n = 0; n < adjList_getNodeNumber(aL); n++) {
        score += getSumOfConsistentAdjacenciesScore(n, aL, ref);
        score += getSumOfConsistentAdjacenciesScore(-n, aL, ref);
    }
    return score;
}
