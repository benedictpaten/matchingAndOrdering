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

refEdge refEdge_construct(int64_t to, double weight) {
    refEdge e;
    e.to = to;
    e.weight = weight;
    return e;
}

int64_t refEdge_to(refEdge *e) {
    return e->to;
}

double refEdge_weight(refEdge *e) {
    return e->weight;
}

static refEdge *refEdge_copy(refEdge *e) {
    refEdge *e2 = st_malloc(sizeof(refEdge));
    *e2 = *e;
    return e2;
}

static int refEdge_cmpByNode(refEdge *e1, refEdge *e2) {
    return e1->to > e2->to ? 1 : (e1->to < e2->to ? -1 : 0);
}

static int refEdge_cmpByWeight(refEdge *e1, refEdge *e2) {
    return e1->weight > e2->weight ? 1 : (e1->weight < e2->weight ? -1 : refEdge_cmpByNode(e1, e2));
}

static int refEdge_cmpByReferencePosition(refEdge *e1, refEdge *e2, reference *ref) {
    return reference_cmp(ref, refEdge_to(e1), refEdge_to(e2));
}

struct _refAdjList {
    stHash **edgeHashes;
    int64_t nodeNumber;
};

refAdjList *refAdjList_construct(int64_t nodeNumber) {
    refAdjList *aL = st_malloc(sizeof(refAdjList));
    aL->nodeNumber = nodeNumber;
    aL->edgeHashes = st_malloc(sizeof(stHash *) * nodeNumber * 2);
    for (int64_t i = 0; i < 2 * nodeNumber; i++) {
        aL->edgeHashes[i] = stHash_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey,
                (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct, free);
    }
    return aL;
}

void refAdjList_destruct(refAdjList *aL) {
    for (int64_t i = 0; i < 2 * aL->nodeNumber; i++) {
        stHash_destruct(aL->edgeHashes[i]);
    }
    free(aL->edgeHashes);
    free(aL);
}

int64_t refAdjList_getNodeNumber(refAdjList *aL) {
    return aL->nodeNumber;
}

static void checkN(int64_t n, int64_t nodeNumber) {
    assert(n <= nodeNumber);
    assert(n != 0);
    assert(n >= -nodeNumber);
}

static int64_t convertN(refAdjList *aL, int64_t n) {
    checkN(n, aL->nodeNumber);
    int64_t i = (n < 0 ? -n : n + aL->nodeNumber) - 1;
    assert(i >= 0);
    assert(i < 2 * aL->nodeNumber);
    return i;
}

double refAdjList_getWeight(refAdjList *aL, int64_t n1, int64_t n2) {
    checkN(n2, aL->nodeNumber);
    stIntTuple *i = stIntTuple_construct1(n2);
    double *weight = stHash_search(aL->edgeHashes[convertN(aL, n1)], i);
    stIntTuple_destruct(i);
    return weight == NULL ? 0.0 : weight[0];
}

static void refAdjList_setWeightP(refAdjList *aL, int64_t n1, int64_t n2, double weight, bool addToWeight) {
    stHash *edges = aL->edgeHashes[convertN(aL, n1)];
    checkN(n2, aL->nodeNumber);
    stIntTuple *i = stIntTuple_construct1(n2);
    double *w = stHash_search(edges, i);
    if (w == NULL) {
        w = st_malloc(sizeof(double));
        stHash_insert(edges, i, w);
        w[0] = weight;
    } else {
        stIntTuple_destruct(i);
        w[0] = addToWeight ? w[0] + weight : weight;
    }
}

void refAdjList_setWeight(refAdjList *aL, int64_t n1, int64_t n2, double weight) {
    refAdjList_setWeightP(aL, n1, n2, weight, 0);
    refAdjList_setWeightP(aL, n2, n1, weight, 0);
}

void refAdjList_addToWeight(refAdjList *aL, int64_t n1, int64_t n2, double weight) {
    refAdjList_setWeightP(aL, n1, n2, weight, 1);
    refAdjList_setWeightP(aL, n2, n1, weight, 1);
}

refAdjListIt adjList_getEdgeIt(refAdjList *aL, int64_t node) {
    refAdjListIt it;
    it.hash = aL->edgeHashes[convertN(aL, node)];
    it.it = stHash_getIterator(it.hash);
    return it;
}

refEdge refAdjListIt_getNext(refAdjListIt *it) {
    refEdge e;
    stIntTuple *i = stHash_getNext(it->it);
    if (i == NULL) {
        e.to = INT64_MAX;
        e.weight = INT64_MAX;
    } else {
        e.to = stIntTuple_get(i, 0);
        e.weight = *(double *) stHash_search(it->hash, i);
    }
    return e;
}

void refAdjListIt_destruct(refAdjListIt *it) {
    stHash_destructIterator(it->it);
}

long double refAdjList_getWeightOfIncidentEdges(refAdjList *aL, int64_t n) {
    long double score = 0.0;
    refAdjListIt it = adjList_getEdgeIt(aL, n);
    refEdge e = refAdjListIt_getNext(&it);
    while (refEdge_to(&e) != INT64_MAX) {
        score += (refEdge_to(&e) == n ? 2 : 1) * refEdge_weight(&e); //Doubles weight of self edges.
        e = refAdjListIt_getNext(&it);
    }
    refAdjListIt_destruct(&it);
    return score;
}

long double refAdjList_getMaxPossibleScore(refAdjList *aL) {
    long double score = 0.0;
    for (int64_t n = 1; n <= refAdjList_getNodeNumber(aL); n++) {
        score += refAdjList_getWeightOfIncidentEdges(aL, n);
        score += refAdjList_getWeightOfIncidentEdges(aL, -n);
    }
    return score / 2.0;
}

int64_t refAdjList_getNumberOfIncidentEdges(refAdjList *aL, int64_t n) {
    int64_t numberOfWeights = 0;
    refAdjListIt it = adjList_getEdgeIt(aL, n);
    refEdge e = refAdjListIt_getNext(&it);
    while (refEdge_to(&e) != INT64_MAX) {
        numberOfWeights++;
        e = refAdjListIt_getNext(&it);
    }
    refAdjListIt_destruct(&it);
    return numberOfWeights;
}

int64_t refAdjList_getNumberOfWeights(refAdjList *aL) {
    int64_t weights = 0;
    for (int64_t n = 1; n <= refAdjList_getNodeNumber(aL); n++) {
        weights += refAdjList_getNumberOfIncidentEdges(aL, n);
        weights += refAdjList_getNumberOfIncidentEdges(aL, -n);
    }
    return weights / 2.0;
}

/*double calculateZScore(int64_t n, int64_t m, int64_t k, double theta) {
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

typedef struct _insertPoint insertPoint;

struct _insertPoint {
    int64_t node;
    int64_t adjNode;
    bool previous;
    long double score;
    bool equivalentInsertPoints;
};

static insertPoint *insertPoint_construct(int64_t node, int64_t adjNode, bool previous, long double score,
        bool equivalentInsertPoints) {
    insertPoint *iP = st_malloc(sizeof(insertPoint));
    iP->node = node;
    iP->adjNode = adjNode;
    iP->score = score;
    iP->previous = previous;
    iP->equivalentInsertPoints = equivalentInsertPoints;
    return iP;
}

static int64_t insertPoint_node(insertPoint *iP) {
    return iP->node;
}

static int64_t insertPoint_adjNode(insertPoint *iP) {
    return iP->adjNode;
}

static long double insertPoint_score(insertPoint *iP) {
    return iP->score;
}

static double insertPoint_previous(insertPoint *iP) {
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
    int64_t node;
    int64_t index;
};

struct _reference {
    int64_t nodeNumber;
    referenceTerm **nodesInGraph;
    stList *referenceIntervals;
};

reference *reference_construct(int64_t nodeNumber) {
    reference *ref = st_malloc(sizeof(reference));
    ref->nodeNumber = nodeNumber;
    ref->nodesInGraph = st_calloc(nodeNumber, sizeof(referenceTerm *));
    ref->referenceIntervals = stList_construct();
    return ref;
}

void reference_destruct(reference *ref) {
    for(int64_t i=0; i<ref->nodeNumber; i++) {
        if(ref->nodesInGraph[i] != NULL) {
            free(ref->nodesInGraph[i]);
        }
    }
    free(ref->nodesInGraph);
    stList_destruct(ref->referenceIntervals);
    free(ref);
}

static referenceTerm *reference_getTerm(reference *ref, int64_t n) {
    assert(llabs(n) <= ref->nodeNumber);
    assert(llabs(n)-1 >= 0);
    return ref->nodesInGraph[llabs(n)-1]; 
}

void reference_insertNodeP(reference *ref, referenceTerm *rT) {
    assert(reference_getTerm(ref, rT->node) == NULL);
    ref->nodesInGraph[llabs(rT->node)-1] = rT;
}

void reference_makeNewInterval(reference *ref, int64_t firstNode, int64_t lastNode) {
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
    reference_insertNodeP(ref, rTF);
    reference_insertNodeP(ref, rTL);
    stList_append(ref->referenceIntervals, rTF);
}

void reference_insertNode(reference *ref, int64_t pNode, int64_t node) {
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
    reference_insertNodeP(ref, rT);
    //Deal with indices
    assert(rT->nTerm->index - rTP->index >= 1);
    if (rT->nTerm->index - rTP->index == 1) { //Need to rebalance
        //Work out the length of the chain
        referenceTerm *rT2 = rTP->first;
        int64_t length = 0;
        while (rT2 != NULL) {
            length++;
            rT2 = rT2->nTerm;
        }
        //Now give every one equi-distant labels.
        assert(length > 1);
        st_logDebug("Rebalancing a reference string with %" PRIi64 " elements\n", length);
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
    reference_insertNode(ref,
            insertPoint_previous(iP) ? insertPoint_adjNode(iP) : reference_getPrevious(ref, insertPoint_adjNode(iP)),
            insertPoint_node(iP));
}

static void reference_removeNode(reference *ref, int64_t n) {
    referenceTerm *rT = reference_getTerm(ref, n);
    assert(rT != NULL);
    if (rT->pTerm == NULL || rT->nTerm == NULL) {
        return;
    }
    ref->nodesInGraph[llabs(n)-1] = NULL;
    rT->nTerm->pTerm = rT->pTerm;
    rT->pTerm->nTerm = rT->nTerm;
    free(rT);
}

bool reference_inGraph(reference *ref, int64_t n) {
    return reference_getTerm(ref, n) != NULL;
}

int64_t reference_getFirstOfInterval(reference *ref, int64_t interval) {
    return ((referenceTerm *) stList_get(ref->referenceIntervals, interval))->node;
}

int64_t reference_getIntervalNumber(reference *ref) {
    return stList_length(ref->referenceIntervals);
}

int64_t reference_getFirst(reference *ref, int64_t n) {
    assert(reference_inGraph(ref, n));
    return reference_getTerm(ref, n)->first->node;
}

int64_t reference_getPrevious(reference *ref, int64_t n) {
    assert(reference_inGraph(ref, n));
    referenceTerm *rT = reference_getTerm(ref, n);
    return rT->pTerm != NULL ? rT->pTerm->node : INT64_MAX;
}

bool reference_getOrientation(reference *ref, int64_t n) {
    assert(reference_inGraph(ref, n));
    return reference_getTerm(ref, n)->node == n;
}

int64_t reference_getNext(reference *ref, int64_t n) {
    assert(reference_inGraph(ref, n));
    referenceTerm *rT = reference_getTerm(ref, n);
    return rT->nTerm != NULL ? rT->nTerm->node : INT64_MAX;
}

int64_t reference_getLast(reference *ref, int64_t n) {
    assert(reference_inGraph(ref, n));
    while (reference_getNext(ref, n) != INT64_MAX) {
        n = reference_getNext(ref, n);
    }
    return n;
}

bool reference_isConsistent(reference *ref, int64_t m, int64_t n) {
    if(reference_getFirst(ref, m) == reference_getFirst(ref,
            n)) {
        if(!reference_getOrientation(ref, m)) {
            return reference_getOrientation(ref, n) && reference_cmp(ref, m, n) == -1;
        }
        return !reference_getOrientation(ref, n) && reference_cmp(ref, n, m) == -1;
    }
    return 0;
}

int reference_cmp(reference *ref, int64_t n1, int64_t n2) {
    referenceTerm *rT1 = reference_getTerm(ref, n1), *rT2 = reference_getTerm(ref, n2);
    assert(rT1 != NULL);
    assert(rT2 != NULL);
    if (rT1->first != rT2->first) {
        return rT1->first > rT2->first ? 1 : -1;
    }
    return rT1->index > rT2->index ? 1 : rT1->index < rT2->index ? -1 : 0;
}

int64_t reference_getRemainingIntervalLength(reference *ref, int64_t n) {
    int64_t j=0;
    while(n != INT64_MAX) {
        j++;
        n = reference_getNext(ref, n);
    }
    return j;
}

void reference_log(reference *ref) {
    st_logInfo("Logging reference with %" PRIi64 " intervals\n", reference_getIntervalNumber(ref));
    for (int64_t i = 0; i < reference_getIntervalNumber(ref); i++) {
        st_logInfo("Interval %" PRIi64 ", nodes:", i);
        int64_t n = reference_getFirstOfInterval(ref, i);
        while (n != INT64_MAX) {
            st_logInfo(" %" PRIi64 ", ", n);
            n = reference_getNext(ref, n);
        }
        st_logInfo("\n");
    }
}

/*
 * Reference algorithm
 */

static stList *getRelevantEdges(refAdjList *aL, reference *ref, int64_t n) {
    stList *edges = stList_construct3(0, free);
    refAdjListIt it = adjList_getEdgeIt(aL, n);
    refEdge e = refAdjListIt_getNext(&it);
    while (refEdge_to(&e) != INT64_MAX) {
        //Check if edge is in graph
        if (reference_inGraph(ref, refEdge_to(&e))) {
            //Is in graph, so add it to list
            stList_append(edges, refEdge_copy(&e));
        }
        e = refAdjListIt_getNext(&it);
    }
    refAdjListIt_destruct(&it);
    //Now do sorting to determine ordering
    stList_sort2(edges, (int(*)(const void *, const void *, const void *)) refEdge_cmpByReferencePosition, ref);
    return edges;
}

static stList *getInsertPointsPrevious(int64_t n, stList *edges, reference *ref) {
    stList *insertPoints = stList_construct();
    long double f = 0.0;
    int64_t intervalName = INT64_MAX;
    for (int64_t i = 0; i < stList_length(edges); i++) {
        refEdge *e = stList_get(edges, i);
        if (reference_getFirst(ref, refEdge_to(e)) != intervalName) { //Reset placement
            f = 0.0;
            intervalName = reference_getFirst(ref, refEdge_to(e));
        }
        if (!reference_getOrientation(ref, refEdge_to(e)) && reference_getNext(ref, refEdge_to(e)) != INT64_MAX) {
            f += refEdge_weight(e);
            assert(reference_getNext(ref, refEdge_to(e)) != INT64_MAX);
            stList_append(insertPoints, insertPoint_construct(n, refEdge_to(e), 1, f, reference_getNext(ref, reference_getNext(ref, refEdge_to(e))) == INT64_MAX ? 0 : 1));
        }
    }
    stList_sort2(insertPoints, (int(*)(const void *, const void *, const void *))insertPoint_cmp, ref);
    return insertPoints;
}

static stList *getInsertPointsNext(int64_t n, stList *edges, reference *ref) {
    stList *insertPoints = stList_construct();
    long double f = 0.0;
    int64_t intervalName = INT64_MAX;
    for (int64_t i = stList_length(edges) - 1; i >= 0; i--) {
        refEdge *e = stList_get(edges, i);
        if (reference_getFirst(ref, refEdge_to(e)) != intervalName) { //Reset placement
            f = 0.0;
            intervalName = reference_getFirst(ref, refEdge_to(e));
        }
        if (reference_getOrientation(ref, refEdge_to(e)) && reference_getPrevious(ref, refEdge_to(e)) != INT64_MAX) {
            f += refEdge_weight(e);
            assert(reference_getPrevious(ref, refEdge_to(e)) != INT64_MAX);
            stList_append(insertPoints, insertPoint_construct(n, refEdge_to(e), 0, f, reference_getPrevious(ref, reference_getPrevious(ref, refEdge_to(e))) == INT64_MAX ? 0 : 1));
        }
    }
    stList_sort2(insertPoints, (int(*)(const void *, const void *, const void *))insertPoint_cmp, ref);
    return insertPoints;
}

static void getInsertionPoints(int64_t n, stList *previousEdges, stList *nextEdges, reference *ref, refAdjList *aL, refAdjList *dAL,
        stList *insertPoints) {
    stList *previousInsertPoints = getInsertPointsPrevious(n, previousEdges, ref);
    stList_appendAll(insertPoints, previousInsertPoints);
    stList *nextInsertPoints = getInsertPointsNext(n, nextEdges, ref);
    stList_appendAll(insertPoints, nextInsertPoints);
    stListIterator *previousInsertPointsIterator = stList_getIterator(previousInsertPoints);
    stListIterator *nextInsertPointsIterator = stList_getIterator(nextInsertPoints);
    insertPoint *iPP = stList_getNext(previousInsertPointsIterator);
    insertPoint *iPN = stList_getNext(nextInsertPointsIterator);
    while (iPP != NULL && iPN != NULL) {
        int64_t i = reference_cmp(ref, insertPoint_adjNode(iPP), insertPoint_adjNode(iPN));
        if (i < 0) {
            while(1) {
                insertPoint *iPPN = stList_getNext(previousInsertPointsIterator);
                if(iPPN != NULL && reference_cmp(ref, insertPoint_adjNode(iPPN), insertPoint_adjNode(iPN)) < 0) {
                    iPP = iPPN;
                    continue;
                }
                if (reference_getFirst(ref, insertPoint_adjNode(iPP)) == reference_getFirst(ref, insertPoint_adjNode(iPN))) {
                    bool equivalentInsertPoints = llabs(reference_getNext(ref, insertPoint_adjNode(iPP))) == llabs(insertPoint_adjNode(iPN)) ? 0 : 1; //I don't think the absolute values are necessary
                    double nWL = refAdjList_getWeight(aL, insertPoint_adjNode(iPP), n); //new weight left
                    double nWR =  refAdjList_getWeight(aL, -n, insertPoint_adjNode(iPN)); //new weight right
                    assert(nWL > 0);
                    assert(nWR > 0);
                    bool left = nWL > nWR;
                    if(equivalentInsertPoints) {
                        assert(reference_getNext(ref, insertPoint_adjNode(iPP)) != INT64_MAX);
                        assert(reference_getPrevious(ref, insertPoint_adjNode(iPN)) != INT64_MAX);
                        double eWL = refAdjList_getWeight(dAL, insertPoint_adjNode(iPP), reference_getNext(ref, insertPoint_adjNode(iPP))); //existing weight left
                        double eWR = refAdjList_getWeight(dAL, -reference_getPrevious(ref, insertPoint_adjNode(iPN)), insertPoint_adjNode(iPN)); //existing weight right
                        if(eWL == 0) {
                            left = eWR > 0 || nWL > nWR;
                        }
                        else {
                            left = !(nWR == 0 || nWR > nWL);
                        }
                    }
                    if (left) {
                        stList_append(
                                insertPoints,
                                insertPoint_construct(n, insertPoint_adjNode(iPP), 1,
                                        insertPoint_score(iPP) + insertPoint_score(iPN), equivalentInsertPoints));
                    } else {
                        stList_append(
                                insertPoints,
                                insertPoint_construct(n, insertPoint_adjNode(iPN), 0,
                                        insertPoint_score(iPP) + insertPoint_score(iPN), equivalentInsertPoints));
                    }
                }
                iPP = iPPN;
                iPN = stList_getNext(nextInsertPointsIterator);
                break;
            }
        } else {
            iPN = stList_getNext(nextInsertPointsIterator);
        }
    }
    stList_destructIterator(previousInsertPointsIterator);
    stList_destructIterator(nextInsertPointsIterator);
    stList_destruct(previousInsertPoints);
    stList_destruct(nextInsertPoints);
}

static insertPoint *getABestInsertNode(int64_t n, refAdjList *aL, refAdjList *dAL, reference *ref) {
    assert(!reference_inGraph(ref, n));
    insertPoint *bestIP = NULL;
    //Get list of edges to nodes already in the reference
    stList *previousEdges = getRelevantEdges(aL, ref, n);
    stList *nextEdges = getRelevantEdges(aL, ref, -n);
    //Now do the hardwork of determining the best insertion point
    stList *insertPoints = stList_construct3(0, free);
    getInsertionPoints(n, previousEdges, nextEdges, ref, aL, dAL, insertPoints);
    getInsertionPoints(-n, nextEdges, previousEdges, ref, aL, dAL, insertPoints);
    //Get best insertion point
    for (int64_t i = 0; i < stList_length(insertPoints); i++) {
        insertPoint *iP = stList_get(insertPoints, i);
        if (bestIP == NULL || insertPoint_score(iP) > insertPoint_score(bestIP) || (insertPoint_score(iP)
                == insertPoint_score(bestIP) && !iP->equivalentInsertPoints)) {
            bestIP = iP;
        }
    }
    if (bestIP != NULL) { //Clone it so we can clean everything else up.
        bestIP = insertPoint_construct(bestIP->node, bestIP->adjNode, bestIP->previous, bestIP->score,
                bestIP->equivalentInsertPoints);
    }
    //Cleanup
    stList_destruct(previousEdges);
    stList_destruct(nextEdges);
    stList_destruct(insertPoints);
    return bestIP;
}

static void insertNode(int64_t n, refAdjList *aL, refAdjList *dAL, reference *ref) {
    if (!reference_inGraph(ref, n)) { //Have a node to add
        insertPoint *bestIP = getABestInsertNode(n, aL, dAL, ref);
        if (bestIP == NULL) { //Make up a location
            st_logDebug("Got a node with no edges linking it into the graph\n");
            reference_insertNode(ref, reference_getFirstOfInterval(ref, 0), n);
        } else {
            reference_insertNode2(ref, bestIP);
            free(bestIP);
        }
    }
}

/*
 * Structure for recording segments connected to segments in the reference but not actually in the reference.
 */

typedef struct _connectedNodeEdge {
    refEdge rE;
    double maxWeight;
    double weightOfEdgesInGraph;
    double inconsistentAdjacencyWeight;
} ConnectedNodeEdge;

struct _connectedNodes {
    stSortedSet *byWeight;
    stSortedSet *byNode;
    refAdjList *aL;
    reference *ref;
    int64_t misses;
};

typedef struct _connectedNodes connectedNodes;

static ConnectedNodeEdge *connectedNodeEdge_construct(refEdge *e, connectedNodes *cN) {
    ConnectedNodeEdge *cNE = st_calloc(1, sizeof(ConnectedNodeEdge));
    cNE->rE = *e;
    cNE->maxWeight = refAdjList_getWeightOfIncidentEdges(cN->aL, refEdge_to(e)) + refAdjList_getWeightOfIncidentEdges(
            cN->aL, -refEdge_to(e));
    return cNE;
}

static long double connectedNodeEdge_calculateWeight(ConnectedNodeEdge *cNE) {
    assert(cNE->inconsistentAdjacencyWeight / cNE->weightOfEdgesInGraph <= 1.0001);
    assert(cNE->weightOfEdgesInGraph / cNE->maxWeight <= 1.0001);
    long double i = (cNE->weightOfEdgesInGraph - cNE->inconsistentAdjacencyWeight) / cNE->maxWeight;
    return i;
}

static void connectedNodes_addNode(connectedNodes *cN, int64_t n) {
    /*
     * Adds nodes not in reference to set of connected nodes.
     */
    refAdjListIt it = adjList_getEdgeIt(cN->aL, n);
    refEdge e = refAdjListIt_getNext(&it);
    while (refEdge_to(&e) != INT64_MAX) {
        if (!reference_inGraph(cN->ref, llabs(refEdge_to(&e)))) {
            e = refEdge_construct(llabs(refEdge_to(&e)), refEdge_weight(&e));
            ConnectedNodeEdge *e2 = stSortedSet_search(cN->byNode, &e);
            if (e2 == NULL) {
                e2 = connectedNodeEdge_construct(&e, cN);
                stSortedSet_insert(cN->byNode, e2);
            } else {
                assert(stSortedSet_search(cN->byWeight, e2) == e2);
                stSortedSet_remove(cN->byWeight, e2);
            }
            e2->weightOfEdgesInGraph += refEdge_weight(&e);
            assert((e2->weightOfEdgesInGraph / e2->maxWeight) <= 1.00001);
            ((refEdge *) e2)->weight = connectedNodeEdge_calculateWeight(e2);
            stSortedSet_insert(cN->byWeight, e2);
        }
        e = refAdjListIt_getNext(&it);
    }
    refAdjListIt_destruct(&it);
}

static connectedNodes *connectedNodes_construct(refAdjList *aL, reference *ref) {
    /*
     * Builds the set of nodes connected to nodes in the reference but not currently in the reference.
     */
    connectedNodes *cN = st_malloc(sizeof(connectedNodes));
    cN->byWeight = stSortedSet_construct3((int(*)(const void *, const void *)) refEdge_cmpByWeight, free);
    cN->byNode = stSortedSet_construct3((int(*)(const void *, const void *)) refEdge_cmpByNode, NULL);
    cN->aL = aL;
    cN->ref = ref;
    cN->misses = 0;
    for (int64_t i = 0; i < reference_getIntervalNumber(ref); i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        connectedNodes_addNode(cN, -n);
        n = reference_getNext(ref, n);
        assert(n != INT64_MAX);
        while (reference_getNext(ref, n) != INT64_MAX) {
            connectedNodes_addNode(cN, n);
            connectedNodes_addNode(cN, -n);
            n = reference_getNext(ref, n);
        }
        connectedNodes_addNode(cN, n);
    }
    return cN;
}

static void connectedNodes_destruct(connectedNodes *cN) {
    stSortedSet_destruct(cN->byNode);
    stSortedSet_destruct(cN->byWeight);
    free(cN);
}

static bool connectedNodes_empty(connectedNodes *cN) {
    return stSortedSet_size(cN->byNode) == 0;
}

static insertPoint *connectedNodes_popBestInsert(connectedNodes *cN, refAdjList *aL, refAdjList *dAL, reference *ref, double wiggle) {
    assert(wiggle <= 1.0);
    int64_t i = 0;
    while (1) {
        ConnectedNodeEdge *cNE = stSortedSet_getLast(cN->byWeight);
        stSortedSet_remove(cN->byWeight, cNE);
        insertPoint *iP = getABestInsertNode(refEdge_to((refEdge *) cNE), aL, dAL, ref);
        assert(iP != NULL);
        cNE->inconsistentAdjacencyWeight = cNE->weightOfEdgesInGraph - insertPoint_score(iP);
        assert(insertPoint_score(iP) / cNE->weightOfEdgesInGraph <= 1.0001);
        ((refEdge *) cNE)->weight = connectedNodeEdge_calculateWeight(cNE); // / (iP->equivalentInsertPoints ? 2 : 1); //Division through by the equivalent insert points number means we try to avoid making totally arbitrary ordering decisions about the partial order.
        if (stSortedSet_size(cN->byWeight) == 0 || refEdge_weight(stSortedSet_getLast(cN->byWeight)) * wiggle
                <= refEdge_weight((refEdge *) cNE) || i++ >= stSortedSet_size(cN->byWeight)) {
            stSortedSet_remove(cN->byNode, cNE);
            return iP;
        }
        stSortedSet_insert(cN->byWeight, cNE);
        cN->misses++;
        free(iP);
    }
}

/*
 * Actual algorithms to make the reference.
 */

void makeReferenceGreedily2(refAdjList *aL, refAdjList *dAL, reference *ref, double wiggle) {
    assert(reference_getIntervalNumber(ref) > 0 || refAdjList_getNodeNumber(aL) == 0);
    connectedNodes *cN = connectedNodes_construct(aL, ref);
    //Iterate over the nodes to check any nodes that are not in the reference
    int64_t i = 0;
    for (int64_t n = 1; n <= refAdjList_getNodeNumber(aL); n++) {
        while (!connectedNodes_empty(cN)) {
            insertPoint *iP = connectedNodes_popBestInsert(cN, aL, dAL, ref, wiggle);
            assert(iP != NULL);
            reference_insertNode2(ref, iP);
            connectedNodes_addNode(cN, insertPoint_node(iP));
            connectedNodes_addNode(cN, -insertPoint_node(iP));
            free(iP);
        }
        if(!reference_inGraph(ref, n)) {
            i++;
            insertNode(n, aL, dAL, ref);
            connectedNodes_addNode(cN, n);
            connectedNodes_addNode(cN, -n);
        }
    }
    st_logDebug("We added %" PRIi64 " unconnected nodes to the reference of %" PRIi64 " nodes\n", i, refAdjList_getNodeNumber(aL));
    st_logDebug("We had %" PRIi64 " connected node misses of %" PRIi64 " nodes\n", cN->misses, refAdjList_getNodeNumber(aL));
    connectedNodes_destruct(cN);
}

void updateReferenceGreedily(refAdjList *aL, refAdjList *dAL, reference *ref, int64_t permutations) {
    for (int64_t i = 0; i < permutations; i++) {
        for (int64_t j = 1; j <= refAdjList_getNodeNumber(aL); j++) {
            int64_t n = st_randomInt(1, refAdjList_getNodeNumber(aL) + 1);
            assert(reference_inGraph(ref, n));
            reference_removeNode(ref, n);
            if (!reference_inGraph(ref, n)) {
                insertNode(n, aL, dAL, ref);
                assert(reference_inGraph(ref, n));
            }
        }
    }
}

static bool nudge(int64_t n, refAdjList *dAL, refAdjList *aL, reference *ref, int64_t maxNudge) {
    //Setup the best insertion spot
    int64_t bestInsert = INT64_MAX;
    int64_t k = reference_getPrevious(ref, n);
    int64_t m = reference_getNext(ref, n);
    assert(k != INT64_MAX && m != INT64_MAX);
    double existingAdjacency1 = refAdjList_getWeight(dAL, -k, n) + refAdjList_getWeight(dAL, -n, m);
    double newAdjacency1 = refAdjList_getWeight(dAL, -k, m);
    double bestScore = existingAdjacency1 - newAdjacency1;

    //Traverse left
    m = k;
    k = reference_getPrevious(ref, m);
    int64_t i = 0;
    while (k != INT64_MAX && refAdjList_getWeight(aL, -m, n) == 0 && i < maxNudge) { //There is a legitimate place to insert, and we are not contradicting any existing weights.
        double existingAdjacency2 = refAdjList_getWeight(dAL, -k, m); //Existing weight at insert point
        double newAdjacency2 = refAdjList_getWeight(dAL, -k, n) + refAdjList_getWeight(dAL, -n, m);
        double newScore = newAdjacency2 - existingAdjacency2;
        if (newScore > bestScore) {
            bestScore = newScore;
            bestInsert = k;
        }
        m = k;
        k = reference_getPrevious(ref, k);
        i++;
    }

    //Traverse right
    k = reference_getNext(ref, n);
    m = reference_getNext(ref, k);
    i = 0;
    while (m != INT64_MAX && refAdjList_getWeight(aL, -n, k) == 0 && i < maxNudge) { //There is a legitimate place to insert, and we are not contradicting any existing weights.
        double existingAdjacency2 = refAdjList_getWeight(dAL, -k, m); //Existing weight at insert point
        double newAdjacency2 = refAdjList_getWeight(dAL, -k, n) + refAdjList_getWeight(dAL, -n, m);
        double newScore = newAdjacency2 - existingAdjacency2;
        if (newScore > bestScore) {
            bestScore = newScore;
            bestInsert = k;
        }
        k = m;
        m = reference_getNext(ref, m);
        i++;
    }

    if (bestInsert != INT64_MAX) {
        reference_removeNode(ref, n);
        reference_insertNode(ref, bestInsert, n);
        return 1;
    }
    return 0;
}

void nudgeGreedily(refAdjList *dAL, refAdjList *aL, reference *ref, int64_t permutations, int64_t maxNudge) {
    for (int64_t i = 0; i < permutations; i++) {
        bool madeNudge = 0;
        for (int64_t j = 0; j < reference_getIntervalNumber(ref); j++) {
            int64_t n = reference_getNext(ref, reference_getFirstOfInterval(ref, j));
            assert(n != INT64_MAX);
            int64_t m;
            while ((m = reference_getNext(ref, n)) != INT64_MAX) {
                madeNudge = madeNudge || nudge(n, dAL, aL, ref, maxNudge);
                n = m;
            }
        }
        if(!madeNudge) {
            st_logDebug("Ran out of nudges after %i iterations\n", i);
            break;
        }
    }
}

static long double getSumOfConsistentAdjacenciesScore(int64_t n, refAdjList *aL, reference *ref) {
    long double score = 0.0;
    refAdjListIt it = adjList_getEdgeIt(aL, n);
    refEdge e = refAdjListIt_getNext(&it);
    while (refEdge_to(&e) != INT64_MAX) {
        if (reference_cmp(ref, n, refEdge_to(&e)) == -1 && reference_getFirst(ref, n) == reference_getFirst(ref,
                refEdge_to(&e)) && !reference_getOrientation(ref, n) && reference_getOrientation(ref, refEdge_to(&e))) {
            score += refEdge_weight(&e);
        }
        e = refAdjListIt_getNext(&it);
    }
    refAdjListIt_destruct(&it);
    return score;
}

long double getReferenceScore(refAdjList *aL, reference *ref) {
    long double score = 0.0;
    for (int64_t n = 1; n <= refAdjList_getNodeNumber(aL); n++) {
        score += getSumOfConsistentAdjacenciesScore(n, aL, ref);
        score += getSumOfConsistentAdjacenciesScore(-n, aL, ref);
    }
    return score;
}

int64_t getBadAdjacencyCount(refAdjList *aL, reference *ref) {
    int64_t badAdjacencies = 0;
    for (int64_t i = 0; i < reference_getIntervalNumber(ref); i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        while (n != INT64_MAX) {
            int64_t m = reference_getNext(ref, n);
            if (m != INT64_MAX && refAdjList_getWeight(aL, -n, m) == 0.0) {
                badAdjacencies++;
            }
            n = m;
        }
    }
    return badAdjacencies;
}

static stList *getValidEdges(int64_t n, refAdjList *aL, reference *ref) {
    /*
     * Get edges from n that are consistent with the reference ordering.
     */
    assert(!reference_getOrientation(ref, n));
    stList *validEdges = stList_construct();
    refAdjListIt it = adjList_getEdgeIt(aL, n);
    refEdge e = refAdjListIt_getNext(&it);
    while (refEdge_to(&e) != INT64_MAX) {
        if(reference_isConsistent(ref, n, refEdge_to(&e))) {
        //if (reference_cmp(ref, n, refEdge_to(&e)) == -1 && reference_getFirst(ref, n) == reference_getFirst(ref,
        //        refEdge_to(&e)) && reference_getOrientation(ref, refEdge_to(&e))) {
            //Is a valid edge
            stList_append(validEdges, refEdge_copy(&e));
        }
        e = refAdjListIt_getNext(&it);
    }
    refAdjListIt_destruct(&it);
    return validEdges;
}

static bool visitP(int64_t n, stSortedSet *visited, stSortedSet *visiting, refAdjList *aL, reference *ref,
        stList *ordering, stList *stack) {
    assert(reference_inGraph(ref, n));
    stIntTuple *i = stIntTuple_construct1(n);
    if (stSortedSet_search(visited, i) == NULL) {
        assert(stSortedSet_search(visiting, i) == NULL); //otherwise we have detected a cycle
        stSortedSet_insert(visiting, i);
        stList *validEdges = getValidEdges(-n, aL, ref); //The minus sign is because we seek edges incident with the righthand side of the segment.
        stList_sort(validEdges, (int(*)(const void *, const void *)) refEdge_cmpByWeight);
        stList_reverse(validEdges); //Traverse edges in reverse order of weight. This should be better, as it will ensure the highest weight adjacency appears in the reference, providing that it can be included in the DFS tree.
        stList_append(stack, i);
        stList_append(stack, validEdges);
        return 1;
    } else {
        assert(stSortedSet_search(visiting, i) != NULL);
        stIntTuple_destruct(i);
        return 0;
    }
}

static void visit(int64_t n, stSortedSet *visited, stSortedSet *visiting, refAdjList *aL, reference *ref,
        stList *ordering, stList *stack) {
    /*
     * Do DFS of nodes using edges that are consistent with the graph.
     */
    if (visitP(n, visited, visiting, aL, ref, ordering, stack)) {
        while (1) {
            stList *edges = stList_peek(stack);
            while (stList_length(edges) > 0) {
                refEdge *e = stList_pop(edges);
                int64_t m = e->to;
                free(e);
                if (visitP(m, visited, visiting, aL, ref, ordering, stack)) {
                    goto top;
                }
            }
            stList_destruct(stList_pop(stack));
            stIntTuple *i = stList_pop(stack);
            stSortedSet_insert(visited, i);
            stList_append(ordering, i);
            if (stList_length(stack) == 0) {
                break;
            }
            top: ;
        }
    }
}

static void reorderReferenceIntervalToAvoidBreakpoints(int64_t startNode, refAdjList *aL, reference *ref) {
    /*
     * Create a topological sort of the nodes in the reference interval, choosing to traverse more highly weighted edges first,
     * with the aim of creating fewer edges that are inconsistent.
     */
    stList *stack = stList_construct();
    stSortedSet *visited = stSortedSet_construct3((int(*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
    stSortedSet *visiting = stSortedSet_construct3((int(*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
    stList *ordering = stList_construct();
    visit(reference_getLast(ref, startNode), visited, visiting, aL, ref, ordering, stack); //Add last node first, as constructed in reverse order.
    visit(startNode, visited, visiting, aL, ref, ordering, stack); //Visit the start node first
    stIntTuple *i = stList_pop(ordering); //Remove the first node, as it must appear first
    int64_t n = startNode;
    while (reference_getNext(ref, n) != INT64_MAX) { //add any other nodes to the ordering that are not on a path from the start node
        visit(n, visited, visiting, aL, ref, ordering, stack);
        n = reference_getNext(ref, n);
    }
    stSortedSet_destruct(visited);
    stSortedSet_destruct(visiting);
    stIntTuple_destruct(i);
    //Now rebuild the reference
    n = reference_getNext(ref, startNode); //Remove the old nodes (this doesn't mess with the first and last nodes).
    while (reference_getNext(ref, n) != INT64_MAX) {
        int64_t m = reference_getNext(ref, n);
        reference_removeNode(ref, n);
        n = m;
    }
    n = startNode; //Now add back the nodes in the new order
    assert(reference_getLast(ref, startNode) != INT64_MAX);
    assert(reference_getLast(ref, startNode) != startNode);
    assert(reference_getNext(ref, startNode) == reference_getLast(ref, startNode));
    i = stList_pop(ordering);
    while (stList_length(ordering) >= 1) {
        int64_t m = stIntTuple_get(i, 0);
        assert(!reference_inGraph(ref, m));
        reference_insertNode(ref, n, m);
        n = m;
        stIntTuple_destruct(i);
        i = stList_pop(ordering);
    }
    stIntTuple_destruct(i);
    stList_destruct(ordering);
}

void reorderReferenceToAvoidBreakpoints(refAdjList *aL, reference *ref) {
    for (int64_t i = 0; i < reference_getIntervalNumber(ref); i++) {
        reorderReferenceIntervalToAvoidBreakpoints(reference_getFirstOfInterval(ref, i), aL, ref);
    }
}

void translocateInterval(reference *ref, int64_t n, int64_t length, int64_t pNode) {
    assert(pNode != INT64_MAX);
    assert(n != INT64_MAX);
    int64_t i=0;
    while(i++<length) {
        int64_t m = reference_getNext(ref, n);
        assert(m != INT64_MAX);
        reference_removeNode(ref, n);
        reference_insertNode(ref, pNode, n);
        pNode = n;
        n = m;
    }
}

int64_t getShortestInterval(reference *ref) {
    int64_t shortestInterval = INT64_MAX;
    int64_t shortestIntervalLength = INT64_MAX;
    for(int64_t i=0; i<reference_getIntervalNumber(ref); i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        int64_t length = reference_getRemainingIntervalLength(ref, n);
        if(length == 0) {
            return n;
        }
        if(length < shortestIntervalLength) {
            shortestInterval = n;
            shortestIntervalLength = length;
        }
    }
    return shortestInterval;
}

void reorderToAvoidOverlargeChromosome(reference *ref, bool (*tooLarge)(reference *, int64_t n)) {
    for(int64_t i=0; i<reference_getIntervalNumber(ref); i++) {
        int64_t n = reference_getFirstOfInterval(ref, i);
        int64_t length = reference_getRemainingIntervalLength(ref, n);
        if(tooLarge(ref, n) && length > 2) {
            int64_t m = getShortestInterval(ref);
            if(m != n) {
                translocateInterval(ref, reference_getNext(ref, n), length/2, m);
            }
        }
    }
}
