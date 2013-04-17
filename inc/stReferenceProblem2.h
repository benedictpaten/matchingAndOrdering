/*
 * referenceProblem.h
 *
 *  Created on: 15 Jan 2013
 *      Author: benedictpaten
 */

#ifndef REFERENCEPROBLEM2_H_
#define REFERENCEPROBLEM2_H_

#include "sonLib.h"

typedef struct _refEdge refEdge;

typedef struct _refAdjList refAdjList;

typedef struct _refAdjListIt refAdjListIt;

typedef struct _reference reference;

typedef struct _referenceTerm referenceTerm;

typedef struct _insertPoint insertPoint;

typedef struct _connectedNodes connectedNodes;

struct _refEdge {
    int64_t to;
    float weight;
};

struct _refAdjList {
    stHash **edgeHashes;
    int64_t nodeNumber;
};

struct _refAdjListIt {
    stHash *hash;
    stHashIterator *it;
};

struct _referenceTerm {
    referenceTerm *first, *pTerm, *nTerm;
    int64_t node;
    int64_t index;
};

struct _reference {
    stHash *nodesInGraph;
    stList *referenceIntervals;
};

struct _insertPoint {
    int64_t node;
    int64_t adjNode;
    bool previous;
    float score;
};

struct _connectedNodes {
    stSortedSet *byWeight;
    stSortedSet *byNode;
    refAdjList *aL;
    reference *ref;
};

/*
 * Adjacency list structure/edge structure
 */
//Negative value indicates the 3' end of segment, positive value indicates 5 prime end.

refEdge edge_construct(int64_t to, float weight);

int64_t edge_to(refEdge *e);

float edge_weight(refEdge *e);

refAdjList *adjList_construct(int64_t nodeNumber);

void adjList_destruct(refAdjList *aL);

int64_t adjList_getNodeNumber(refAdjList *aL);

float adjList_getWeight(refAdjList *aL, int64_t n1, int64_t n2);

void adjList_setWeight(refAdjList *aL, int64_t n1, int64_t n2, float weight);

void adjList_addToWeight(refAdjList *aL, int64_t n1, int64_t n2, float weight);

refAdjListIt adjList_getEdgeIt(refAdjList *aL, int64_t node);

refEdge adjListIt_getNext(refAdjListIt *it);

void adjListIt_destruct(refAdjListIt *it);

double adjList_getMaxPossibleScore(refAdjList *aL);

//double calculateZScore(int64_t n, int64_t m, int64_t k, double theta);

/*
 * Reference structure
 */

reference *reference_construct();

void reference_destruct(reference *ref);

//Need to make one or more intervals before you can insert other nodes into reference.
void reference_makeNewInterval(reference *ref, int64_t leftNode, int64_t rightNode);

void reference_insertNode(reference *ref, int64_t pNode, int64_t node);

bool reference_inGraph(reference *ref, int64_t n);

int64_t reference_getFirstOfInterval(reference *ref, int64_t interval);

int64_t reference_getIntervalNumber(reference *ref);

int64_t reference_getFirst(reference *ref, int64_t n);

int64_t reference_getPrevious(reference *ref, int64_t n);

int64_t reference_getLast(reference *ref, int64_t n);

//Returns nonzero if segment is in reference in same orientation, traversing the reference sequence(s) from 5' to 3'.
bool reference_getOrientation(reference *ref, int64_t n);

//Gets the next position within the reference.
int64_t reference_getNext(reference *ref, int64_t n);

//Compares segments position within a reference, ignoring orientation.
int reference_cmp(reference *ref, int64_t n1, int64_t n2);

void reference_log(reference *ref);

/*
 * Reference algorithms
 */

void makeReferenceGreedily2(refAdjList *aL, reference *ref);

void updateReferenceGreedily(refAdjList *aL, reference *ref, int64_t permutations);

/*
 * Create a topological sort of each reference interval, trying to place nodes that are connected by direct adjacencies next to one another.
 */
void reorderReferenceToAvoidBreakpoints(refAdjList *aL, reference *ref);

double getReferenceScore(refAdjList *aL, reference *ref);

/*
 * Nudge the blocks to try to eliminate "bad adjacencies", where to blocks are adjacent but do not have a direct weight between them.
 */
void nudgeGreedily(refAdjList *dAL, refAdjList *aL, reference *ref, int64_t permutations, int64_t maxNudge);

/*
 * Count of adjacent nodes in reference that have no edge connecting them.
 */
int64_t getBadAdjacencyCount(refAdjList *aL, reference *ref);

#endif /* REFERENCEPROBLEM2_H_ */
