/*
 * DBGraph_reduced.h
 *
 *  Created on: Nov 14, 2014
 *      Author: kaempfpp
 */

#ifndef DBGRAPH_REDUCED_H_
#define DBGRAPH_REDUCED_H_

int verticalReduction();

void reduceRedGraph();

void collapseEdges(int i, int dest,int ori, int up);

void reduceRedGraph_strong();

int scaffoldRedGraph();

void redGraphConnector();

void printRedGraphToFile(char* filePath);

void printRedGraph();

void printRedGraphEdge();

void printRedGraphList(char*);

void printRedDot();

void travToRedOVL();

void travToRedOVL_v2();

void collectInternalNodes(int i,int Knoten);

int outDegree(int node);

struct Knoten* addKnoten();

void prepareGraph();

void labelGraph();

#endif /* DBGRAPH_REDUCED_H_ */
