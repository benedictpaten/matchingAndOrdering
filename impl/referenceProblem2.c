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
    assert(n <= nodeNumber);
    assert(n != 0);
    assert(n >= -nodeNumber);
}

static int32_t convertN(adjList *aL, int32_t n) {
    checkN(n, aL->nodeNumber);
    int32_t i = (n < 0 ? -n : n + aL->nodeNumber) - 1;
    assert(i >= 0);
    assert(i < 2 * aL->nodeNumber);
    return i;
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
        w = st_malloc(sizeof(float));
        stHash_insert(edges, i, w);
    } else {
        stIntTuple_destruct(i);
    }
    w[0] = weight;
}

void adjList_setWeight(adjList *aL, int32_t n1, int32_t n2, float weight) {
    adjList_setWeightP(aL, n1, n2, weight);
    adjList_setWeightP(aL, n2, n1, weight);
}

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
        if (edge_to(&e) >= n) {
            score += edge_weight(&e);
        }
        e = adjListIt_getNext(&it);
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
    return score;
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
    bool previous;
    float score;
};

typedef struct _insertPoint insertPoint;

static insertPoint *insertPoint_construct(int32_t node, int32_t adjNode, bool previous, float score) {
    insertPoint *insertPoint = st_malloc(sizeof(insertPoint));
    insertPoint->node = node;
    insertPoint->adjNode = adjNode;
    insertPoint->score = score;
    insertPoint->previous = previous;
    return insertPoint;
}

static int32_t insertPoint_node(insertPoint *iP) {
    return iP->node;
}

static int32_t insertPoint_adjNode(insertPoint *iP) {
    return iP->adjNode;
}

static float insertPoint_score(insertPoint *iP) {
    return iP->score;
}

static float insertPoint_previous(insertPoint *iP) {
    return iP->previous;
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

static uint32_t referenceTermHashKey(const void *refTerm) {
    return abs(((referenceTerm *) refTerm)->node);
}

static int referenceTermEqualsFn(const void *refTerm1, const void *refTerm2) {
    return abs(((referenceTerm *) refTerm1)->node) == abs(((referenceTerm *) refTerm2)->node);
}

reference *reference_construct() {
    reference *ref = st_malloc(sizeof(reference));
    ref->nodesInGraph = stHash_construct3(referenceTermHashKey, referenceTermEqualsFn, free, NULL);
    ref->referenceIntervals = stList_construct();
    return ref;
}

void reference_destruct(reference *ref) {
    stHash_destruct(ref->nodesInGraph);
    stList_destruct(ref->referenceIntervals);
    free(ref);
}

void reference_makeNewInterval(reference *ref, int32_t firstNode, int32_t lastNode) {
    referenceTerm *rTF = st_malloc(sizeof(referenceTerm)), *rTL = st_malloc(sizeof(referenceTerm));
    rTF->node = firstNode;
    rTL->node = lastNode;
    rTF->nTerm = rTL;
    rTL->pTerm = rTF;
    rTF->pTerm = NULL;
    rTL->nTerm = NULL;
    rTF->first = rTF;
    rTL->first = rTF;
    rTF->index = 0;
    rTL->index = INT64_MAX; //This forces the rebalancing code to be exercised.
    stHash_insert(ref->nodesInGraph, rTF, rTF);
    stHash_insert(ref->nodesInGraph, rTL, rTL);
    stList_append(ref->referenceIntervals, rTF);
}

static referenceTerm *reference_getTerm(reference *ref, int32_t n) {
    referenceTerm rT;
    rT.node = n;
    return stHash_search(ref->nodesInGraph, &rT);
}

void reference_insertNode(reference *ref, int32_t pNode, int32_t node) {
    referenceTerm *rT = st_malloc(sizeof(referenceTerm)), *rTP;
    rT->node = node;
    rTP = reference_getTerm(ref, pNode);
    assert(rTP != NULL);
    rT->nTerm = rTP->nTerm;
    assert(rT->nTerm != NULL);
    rT->pTerm = rTP;
    rTP->nTerm = rT;
    rT->nTerm->pTerm = rT;
    rT->first = rTP->first;
    stHash_insert(ref->nodesInGraph, rT, rT);
    //Deal with indices
    assert(rT->nTerm->index - rTP->index >= 1);
    if (rT->nTerm->index - rTP->index == 1) { //Need to rebalance
        //Work out the length of the chain
        referenceTerm *rT2 = rTP->first;
        int32_t length = 0;
        while (rT2 != NULL) {
            length++;
            rT2 = rT2->nTerm;
        }
        //Now give every one equi-distant labels.
        assert(length > 1);
        st_logDebug("Rebalancing a reference string with %i elements\n", length);
        int64_t spacer = INT64_MAX / (length - 1);
        assert(spacer > 0);
        rT2 = rTP->first;
        rT2->index = 0;
        while (rT2->nTerm != NULL) {
            rT2->nTerm->index = rT2->index + spacer;
            rT2 = rT2->nTerm;
        }
    }
    assert(rT->nTerm->index - rTP->index > 1);
    rT->index = rTP->index + (rT->nTerm->index - rTP->index) / 2;
}

static void reference_insertNode2(reference *ref, insertPoint *iP) {
    reference_insertNode(ref, insertPoint_previous(iP) ? insertPoint_adjNode(iP) : reference_getPrevious(ref, insertPoint_adjNode(iP)), insertPoint_node(iP));
}

static void reference_removeNode(reference *ref, int32_t n) {
    referenceTerm *rT = reference_getTerm(ref, n);
    assert(rT != NULL);
    if (rT->pTerm == NULL || rT->nTerm == NULL) {
        return;
    }
    stHash_remove(ref->nodesInGraph, rT);
    rT->nTerm->pTerm = rT->pTerm;
    rT->pTerm->nTerm = rT->nTerm;
    free(rT);
}

bool reference_inGraph(reference *ref, int32_t n) {
    return reference_getTerm(ref, n) != NULL;
}

int32_t reference_getFirstOfInterval(reference *ref, int32_t interval) {
    return ((referenceTerm *) stList_get(ref->referenceIntervals, interval))->node;
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
    for (int32_t i = 0; i < reference_getIntervalNumber(ref); i++) {
        st_logInfo("Interval %i, nodes:", i);
        int32_t n = reference_getFirstOfInterval(ref, i);
        while (n != INT32_MAX) {
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
    stList_sort2(edges, (int(*)(const void *, const void *, const void *)) edge_cmp, ref);
    return edges;
}

static void getInsertPointsPrevious(int32_t n, stList *edges, reference *ref, stList *insertPoints) {
    double f = 0.0;
    int32_t intervalName = INT32_MAX;
    for (int32_t i = 0; i < stList_length(edges); i++) {
        edge *e = stList_get(edges, i);
        if (reference_getFirst(ref, edge_to(e)) != intervalName) { //Reset placement
            f = 0.0;
            intervalName = reference_getFirst(ref, edge_to(e));
        }
        if (!reference_getOrientation(ref, edge_to(e)) && reference_getNext(ref, edge_to(e)) != INT32_MAX) {
            f += edge_weight(e);
            stList_append(insertPoints, insertPoint_construct(n, edge_to(e), 1, f));
        }
    }
}

static void getInsertPointsNext(int32_t n, stList *edges, reference *ref, stList *insertPoints) {
    double f = 0.0;
    int32_t intervalName = INT32_MAX;
    for (int32_t i = stList_length(edges) - 1; i >= 0; i--) {
        edge *e = stList_get(edges, i);
        if (reference_getFirst(ref, edge_to(e)) != intervalName) { //Reset placement
            f = 0.0;
            intervalName = reference_getFirst(ref, edge_to(e));
        }
        if (reference_getOrientation(ref, edge_to(e)) && reference_getPrevious(ref, edge_to(e)) != INT32_MAX) {
            f += edge_weight(e);
            stList_append(insertPoints, insertPoint_construct(n, edge_to(e), 0, f));
        }
    }
}

static void getInsertionPoints(int32_t n, stList *previousEdges, stList *nextEdges, reference *ref, adjList *aL, stList *insertPoints) {
    stList *insertPoints2 = stList_construct();
    getInsertPointsPrevious(n, previousEdges, ref, insertPoints2);
    getInsertPointsNext(n, nextEdges, ref, insertPoints2);
    stList_sort2(insertPoints2, (int(*)(const void *, const void *, const void *)) insertPoint_cmp, ref);
    for (int32_t i = 1; i < stList_length(insertPoints2); i++) {
        insertPoint *iPP = stList_get(insertPoints2, i - 1);
        insertPoint *iPN = stList_get(insertPoints2, i);
        if (reference_getFirst(ref, insertPoint_adjNode(iPP))
                == reference_getFirst(ref, insertPoint_adjNode(iPN)) &&
                insertPoint_previous(iPP) && !insertPoint_previous(iPN)) {
            if (adjList_getWeight(aL, n, insertPoint_adjNode(iPP)) > adjList_getWeight(aL, -n, insertPoint_adjNode(iPN))) {
                stList_append(insertPoints,
                        insertPoint_construct(n, insertPoint_adjNode(iPP), 1, insertPoint_score(iPP) + insertPoint_score(iPN)));
            } else {
                stList_append(insertPoints,
                        insertPoint_construct(n, insertPoint_adjNode(iPN), 0, insertPoint_score(iPP) + insertPoint_score(iPN)));
            }
        }
    }
    stList_appendAll(insertPoints, insertPoints2);
    stList_destruct(insertPoints2);
}

static void insertNode(int32_t n, adjList *aL, reference *ref) {
    if (!reference_inGraph(ref, n)) { //Have a node to add
        //Get list of edges to nodes already in the reference
        stList *previousEdges = getRelevantEdges(aL, ref, n);
        stList *nextEdges = getRelevantEdges(aL, ref, -n);
        //Now do the hardwork of determining the best insertion point
        stList *insertPoints = stList_construct3(0, free);
        getInsertionPoints(n, previousEdges, nextEdges, ref, aL, insertPoints);
        getInsertionPoints(-n, nextEdges, previousEdges, ref, aL, insertPoints);
        //Get best insertion point
        insertPoint *bestIP = NULL;
        for (int32_t i = 0; i < stList_length(insertPoints); i++) {
            insertPoint *iP = stList_get(insertPoints, i);
            if (bestIP == NULL || insertPoint_score(iP) > insertPoint_score(bestIP)) {
                bestIP = iP;
            }
        }
        if(bestIP == NULL) { //Make up a location
            st_logDebug("Got a node with no edges linking it into the graph\n");
            reference_insertNode(ref, reference_getFirstOfInterval(ref, 0), n);
        }
        else {
            reference_insertNode2(ref, bestIP);
        }
        //Cleanup
        stList_destruct(previousEdges);
        stList_destruct(nextEdges);
        stList_destruct(insertPoints);
    }
}

void makeReferenceGreedily2(adjList *aL, reference *ref) {
    assert(reference_getIntervalNumber(ref) > 0 || adjList_getNodeNumber(aL) == 0);
    for (int32_t n = 1; n <= adjList_getNodeNumber(aL); n++) {
        insertNode(n, aL, ref);
    }
}

void updateReferenceGreedily(adjList *aL, reference *ref, int32_t permutations) {
    for (int32_t i = 0; i < permutations; i++) {
        for (int32_t j = 1; j <= adjList_getNodeNumber(aL); j++) {
            int32_t n = st_randomInt(1, adjList_getNodeNumber(aL)+1);
            assert(reference_inGraph(ref, n));
            reference_removeNode(ref, n);
            if (!reference_inGraph(ref, n)) {
                insertNode(n, aL, ref);
                assert(reference_inGraph(ref, n));
            }
        }
    }
}

static double getSumOfConsistentAdjacenciesScore(int32_t n, adjList *aL, reference *ref) {
    double score = 0.0;
    adjListIt it = adjList_getEdgeIt(aL, n);
    edge e = adjListIt_getNext(&it);
    while (edge_to(&e) != INT32_MAX) {
        if (reference_cmp(ref, n, edge_to(&e)) == -1 &&
                reference_getFirst(ref, n) == reference_getFirst(ref, edge_to(&e)) &&
                !reference_getOrientation(ref, n) && reference_getOrientation(ref, edge_to(&e))) {
            score += edge_weight(&e);
        }
        e = adjListIt_getNext(&it);
    }
    adjListIt_destruct(&it);
    return score;
}

double getReferenceScore(adjList *aL, reference *ref) {
    double score = 0.0;
    for (int32_t n = 1; n <= adjList_getNodeNumber(aL); n++) {
        score += getSumOfConsistentAdjacenciesScore(n, aL, ref);
        score += getSumOfConsistentAdjacenciesScore(-n, aL, ref);
    }
    return score;
}