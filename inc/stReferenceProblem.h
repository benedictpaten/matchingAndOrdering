/*
 * adjacencyProblem.h
 *
 *  Created on: 11 Aug 2011
 *      Author: benedictpaten
 */

#ifndef REFERENCEPROBLEM_H_
#define REFERENCEPROBLEM_H_

stList *makeReferenceGreedily(stList *stubs, stList *chains,
        float *z, int64_t nodeNumber, double *totalScore, bool fast);

void greedyImprovement(stList *reference, stList *chains, float *z, int64_t permutations);

void gibbsSamplingWithSimulatedAnnealing(stList *reference,
        stList *chains, float *z, int64_t permutations,
        double(*temperature)(double));

stList *convertReferenceToAdjacencyEdges(stList *reference);

double exponentiallyDecreasingTemperatureFn(double d);

double constantTemperatureFn(double d);

void logReference(stList *reference, int64_t nodeNumber, float *zMatrix, double totalScore, const char *message);

double calculateZScoreOfReference(stList *reference, int64_t nodeNumber, float *zMatrix);

void logZScoreOfReference(stList *reference, int64_t nodeNumber, float *zMatrix);

double calculateMaxZ(int64_t nodeNumber, float *zMatrix);

long double calculateZScore(int64_t n, int64_t m, int64_t k, long double theta);

#endif /* REFERENCEPROBLEM_H_ */
