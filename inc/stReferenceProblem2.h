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

struct _refEdge {
    int64_t to;
    double weight;
};

struct _refAdjListIt {
    stHash *hash;
    stHashIterator *it;
};

/*
 * Adjacency list structure/edge structure
 */
//Negative value indicates the 3' end of segment, positive value indicates 5 prime end.

refEdge refEdge_construct(int64_t to, double weight);

int64_t refEdge_to(refEdge *e);

double refEdge_weight(refEdge *e);

refAdjList *refAdjList_construct(int64_t nodeNumber);

void refAdjList_destruct(refAdjList *aL);

int64_t refAdjList_getNodeNumber(refAdjList *aL);

double refAdjList_getWeight(refAdjList *aL, int64_t n1, int64_t n2);

void refAdjList_setWeight(refAdjList *aL, int64_t n1, int64_t n2, double weight);

void refAdjList_addToWeight(refAdjList *aL, int64_t n1, int64_t n2, double weight);

refAdjListIt adjList_getEdgeIt(refAdjList *aL, int64_t node);

refEdge refAdjListIt_getNext(refAdjListIt *it);

void refAdjListIt_destruct(refAdjListIt *it);

double refAdjList_getMaxPossibleScore(refAdjList *aL);

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
