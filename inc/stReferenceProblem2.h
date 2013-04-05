/*
 * referenceProblem.h
 *
 *  Created on: 15 Jan 2013
 *      Author: benedictpaten
 */

#ifndef REFERENCEPROBLEM2_H_
#define REFERENCEPROBLEM2_H_

#include "sonLib.h"

struct _edge {
    int32_t to;
    float weight;
};

typedef struct _edge edge;

typedef struct _adjList adjList;

struct _adjListIt {
    stHash *hash;
    stHashIterator *it;
};

typedef struct _adjListIt adjListIt;

typedef struct _reference reference;

/*
 * Adjacency list structure/edge structure
 */
//Negative value indicates the 3' end of segment, positive value indicates 5 prime end.

edge edge_construct(int32_t to, float weight);

int32_t edge_to(edge *e);

float edge_weight(edge *e);

adjList *adjList_construct(int32_t nodeNumber);

void adjList_destruct(adjList *aL);

int32_t adjList_getNodeNumber(adjList *aL);

float adjList_getWeight(adjList *aL, int32_t n1, int32_t n2);

void adjList_setWeight(adjList *aL, int32_t n1, int32_t n2, float weight);

void adjList_addToWeight(adjList *aL, int32_t n1, int32_t n2, float weight);

adjListIt adjList_getEdgeIt(adjList *aL, int32_t node);

edge adjListIt_getNext(adjListIt *it);

void adjListIt_destruct(adjListIt *it);

double adjList_getMaxPossibleScore(adjList *aL);

//double calculateZScore(int32_t n, int32_t m, int32_t k, double theta);

/*
 * Reference structure
 */

reference *reference_construct();

void reference_destruct(reference *ref);

//Need to make one or more intervals before you can insert other nodes into reference.
void reference_makeNewInterval(reference *ref, int32_t leftNode, int32_t rightNode);

void reference_insertNode(reference *ref, int32_t pNode, int32_t node);

bool reference_inGraph(reference *ref, int32_t n);

int32_t reference_getFirstOfInterval(reference *ref, int32_t interval);

int32_t reference_getIntervalNumber(reference *ref);

int32_t reference_getFirst(reference *ref, int32_t n);

int32_t reference_getPrevious(reference *ref, int32_t n);

int32_t reference_getLast(reference *ref, int32_t n);

//Returns nonzero if segment is in reference in same orientation, traversing the reference sequence(s) from 5' to 3'.
bool reference_getOrientation(reference *ref, int32_t n);

//Gets the next position within the reference.
int32_t reference_getNext(reference *ref, int32_t n);

//Compares segments position within a reference, ignoring orientation.
int reference_cmp(reference *ref, int32_t n1, int32_t n2);

void reference_log(reference *ref);

/*
 * Reference algorithms
 */

void makeReferenceGreedily2(adjList *aL, reference *ref);

void updateReferenceGreedily(adjList *aL, reference *ref, int32_t permutations);

/*
 * Create a topological sort of each reference interval, trying to place nodes that are connected by direct adjacencies next to one another.
 */
void reorderReferenceToAvoidBreakpoints(adjList *aL, reference *ref);

double getReferenceScore(adjList *aL, reference *ref);

#endif /* REFERENCEPROBLEM2_H_ */
