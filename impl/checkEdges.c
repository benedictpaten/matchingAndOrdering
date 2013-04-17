#include "sonLib.h"
#include "shared.h"

void checkEdges(stList *edges, stSortedSet * nodes, bool coversAllNodes,
        bool isClique) {
    int64_t nodeNumber = stSortedSet_size(nodes);
    int64_t maxNode = nodeNumber == 0 ? 0 : stIntTuple_get(stSortedSet_getLast(nodes), 0);
    (void)maxNode;
    stSortedSet *edgesSeen = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        /*
         * Check edges connect actual nodes.
         */
        assert(stIntTuple_get(edge, 0) >= 0);
        assert(stIntTuple_get(edge, 0) <= maxNode);
        assert(stIntTuple_get(edge, 1) >= 0);
        assert(stIntTuple_get(edge, 1) <= maxNode);
        assert(
                stIntTuple_get(edge, 1) != stIntTuple_get(edge,
                        0)); //No self edges!
        /*
         * Check is not a multi-graph.
         */
        assert(
                !edgeInSet(edgesSeen, stIntTuple_get(edge, 0),
                        stIntTuple_get(edge, 1)));
        addEdgeToSet(edgesSeen, stIntTuple_get(edge, 0),
                stIntTuple_get(edge, 1));
        /*
         * Check weight, if  weighted.
         */
        if (stIntTuple_length(edge) == 3) {
            assert(stIntTuple_get(edge, 2) >= 0);
        } else {
            assert(stIntTuple_length(edge) == 2);
        }
    }
    stSortedSet_destruct(edgesSeen);
    /*
     * Check edges connect all nodes.
     */
    if (coversAllNodes) {
        stSortedSet *nodes = getNodeSetOfEdges(edges);
        assert(stSortedSet_size(nodes) == nodeNumber);
        stSortedSet_destruct(nodes);
    }
    /*
     * Check is clique.
     */
    if (isClique) {
        assert(
                stList_length(edges) == (nodeNumber * nodeNumber - nodeNumber)
                        / 2);
    }
}
