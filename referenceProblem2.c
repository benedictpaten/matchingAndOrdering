
/*
 * Simple set of lists
 */
void stAdjacencyList_addEdge(stAdjacencyList *adjacencyList, void *node1, void *node2);

stList *stAdjacencyList_getEdges(stAdjacencyList *adjacencyList, void *node1);

/*
 * Based upon hash
 */

void stSparseMatrix_set(stSparseMatrix *sparseMatrix, void *x, void *y, void *value);

void *stSparseMatrix_get(stSparseMatrix *sparseMatrix, void *x, void *y);


/*
 * Based upon random binary tree.
 */

stBinaryTree *makeChain(void *node1, void *node2);

void stBinaryTree_insertInInterval(stBinaryTree *start, stBinaryTree *end);

int stBinaryTree_compareNodes(void *node1, void *node2);

void stBinaryTree_removeNode(stBinaryTree *node);

void stBinaryTree_setLeftChild(stBinaryTree *parent, stBinaryTree *child);

void stBinaryTree_setRightChild(stBinaryTree *parent, stBinaryTree *child);





void *addNode(stRandomBinaryTree *permutation, void *predecessor, void *object);

void removeNode(stRandomBinaryTree *permutation, void *node);



/*
 *  Create adjacency list / sparse matrix
 *
 *  Create random binary matrix
 *
 *  Create empty set of intervals in binary tree
 *
 *  In random order, for each chain:
 *
 *      Get other nodes connected to chain
 *
 *      Using binary tree, order those nodes
 *
 *      Find maximal insertion point
 *
 *      Add to tree
 */
