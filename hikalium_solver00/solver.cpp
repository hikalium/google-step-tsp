#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#define MAX_NODE_SIZE	(1 << 22)

class Vector2D
{
public:
	double x, y;
	void set(double x, double y)
	{
		this->x = x;
		this->y = y;
	}
	double calcDistanceTo(Vector2D *another)
	{
		double dx = (another->x - x);
		double dy = (another->y - y);
		return sqrt(dx * dx + dy * dy);
	}
	void print()
	{	
		std::cout << x << ", " << y << std::endl;
	}
	static double calcCrossProduct(Vector2D *vL, Vector2D *vR)
	{
		return (vL->x * vR->y) - (vL->y * vR ->x);
	}
	static double calcDotProduct(Vector2D *vL, Vector2D *vR)
	{
		return (vL->x * vR->x) + (vL->y * vR ->y);
	}
	double calcLength()
	{
		return sqrt(x * x + y * y);
	}
	Vector2D getVectorTo(Vector2D *to)
	{
		Vector2D v;
		v.x = to->x - x;
		v.y = to->y - y;
		return v;
	}
	double calcDistanceFromLineSegmentPQ(Vector2D *p, Vector2D *q)
	{
		Vector2D pq = p->getVectorTo(q);
		Vector2D pt = p->getVectorTo(this);
		double s = fabs(calcCrossProduct(&pq, &pt));
		double l = pq.calcLength();
		if(l == 0){
			return getVectorTo(p).calcLength();
		}
		double d = s / l;

		s = calcDotProduct(&pq, &pt); 
		if(s < 0){
			l = - (s / l);
			d = sqrt(d * d + l * l);
		} else if(s > l * l){
			l = s / l - l;
			d = sqrt(d * d + l * l);
		}
		return d;
	}
};

class Node : public Vector2D
{
	typedef Vector2D base;
public:
	int index/*, path_index*/;
	/*double distance_to_path;
	bool used;*/
	Node(){
		/*used = false;*/
	}
	void print()
	{	
		std::cout << index << "\t";
		base::print();
	}
};

class NodePair
{
public:
	Node *p, *q;
	NodePair()
	{
		NodePair(nullptr, nullptr);
	}
	NodePair(Node *p, Node *q)
	{
		this->p = p;
		this->q = q;
	}
	double calcDistance()
	{
		return p->calcDistanceTo(q);
	}
	double calcDistanceFromThisLineSegment(Node *point)
	{
		return point->calcDistanceFromLineSegmentPQ(p, q);
	}
	void print()
	{	
		std::cout << "Pair: " << std::endl;
		p->print();
		q->print();
		std::cout << "distance: " << calcDistance() << std::endl;
	}
};

class NodeList
{
private:
	int size;
	int max_size;
public:
	Node *list;
	NodeList(int max)
	{
		list = new Node[max];
		assert(list);
		max_size = max;
		size = 0;
	}
	int getSize()
	{
		return size;
	}
	Node *getNode(int index)
	{
		assert(0 <= index && index < size);
		return &list[index];
	}
	int append(double x, double y)
	{
		if(size >= max_size) return 1;
		list[size].x = x;
		list[size].y = y;
		list[size].index = size;
		size++;
		return 0;
	}
	int append(Node *node)
	{
		if(size >= max_size) return 1;
		list[size].x = node->x;
		list[size].y = node->y;
		list[size].index = node->index;
		size++;
		return 0;
	}
	void printAll()
	{
		for(int i = 0; i < size; i++){
			list[i].print();
		}
	}
	NodePair getMostDistantPair()
	{
		NodePair maxPair;
		double maxDistance = 0;
		for(int i = 0; i < size; i++){
			for(int k = 0; k < size; k++){
				double distance = list[i].calcDistanceTo(&list[k]);
				if(distance > maxDistance){
					maxDistance = distance;
					maxPair.p = &list[i];
					maxPair.q = &list[k];
				}
			}
		}
		return maxPair;
	}
	int readFromFile(const char *filename)
	{
		FILE *fp;
		fp = fopen(filename, "rb");
		if(!fp){
			std::cout << "input file not found." << std::endl;
			return 0;
		}
		char buf[128];
		fgets(buf, sizeof(buf), fp);	// skip first line
		double x, y;
		try{
			while(!feof(fp)){
				if(fscanf(fp, "%lf,%lf", &x, &y) != 2) break;
				append(x, y);
			}
		} catch(std::invalid_argument){
		}
		// printAll();
		std::cout << "read " << size << " nodes" << std::endl;
		return size;
	}
};

class NodeIndexList
{
private:
	NodeList *base_list;
	int *list;
	int size;
public:
	int max_size;
	NodeIndexList(NodeList *base_list, int size)
	{
		assert(base_list);
		this->base_list = base_list;
		max_size = size;
		list = new int[max_size];
		max_size = size;
		this->size = 0;
	}
	NodeIndexList(NodeIndexList *base_list, int size)
	{
		assert(base_list && base_list->base_list);
		this->base_list = base_list->base_list;
		max_size = size;
		list = new int[max_size];
		this->size = 0;
	}
	int getSize()
	{
		return size;
	}
	Node *getNode(int index)
	{
		return this->base_list->getNode(list[index]);
	}
	NodeList *getBaseList()
	{
		return base_list;
	}
	int append(int index)
	{
		assert(size < max_size);
		list[size] = index;
		size++;
		return 0;
	}
	int pop()
	{
		assert(size > 0);
		size--;
		int index = list[size];
		list[size] = 0;
		return index;
	}
	void appendAll(NodeList *base_list)
	{
		assert(base_list->getSize() <= max_size);
		this->base_list = base_list;
		size = base_list->getSize();
		for(int i = 0; i < size; i++){
			list[i] = base_list->getNode(i)->index;	
		}
	}
	void appendAll(NodeIndexList *base_list)
	{
		assert(base_list->getSize() <= max_size);
		this->base_list = base_list->base_list;
		size = base_list->getSize();
		for(int i = 0; i < size; i++){
			list[i] = base_list->list[i];	
		}
	}
	void printAll()
	{
		for(int i = 0; i < size; i++){
			base_list->getNode(list[i])->print();
		}
	}
	void dumpAllIndex(std::ofstream *ofs)
	{
		for(int i = 0; i < size; i++){
			*ofs << list[i] << std::endl;
			std::cout << list[i] << std::endl;
		}
	}
	bool hasIndex(int index)
	{
		for(int i = 0; i < size; i++){
			if(list[i] == index) return true;
		}
		return false;
	}
	double getLastTwoDistance()
	{
		assert(size >= 2);
		NodePair pair(getNode(size - 1), getNode(size - 2));
		return pair.calcDistance();
	}
	double getDistanceBetweenFirstAndLast()
	{
		assert(size >= 2);
		NodePair pair(getNode(size - 1), getNode(0));
		return pair.calcDistance();
	}
	void saveSolution(std::string filename)
	{
		std::ofstream ofs;
		ofs.open(filename, std::ios::out);
		ofs << "index" << std::endl;
		std::cout << "index" << std::endl;
		for(int i = 0; i < size; i++){
			ofs << list[i] << std::endl;
			//std::cout << list[i]->index << std::endl;
		}
		std::cout << "solution saved to " << filename << std::endl;
	}
};

class BruteForceSolver
{
	// 16 -> 61s
	// 8 -> 0.1s 
	public:
	NodeIndexList *base_list;
	double bestPathLen;
	NodeIndexList *best_list;
	//static const int updateLimit = 10;
	//int updateCount;
	BruteForceSolver(NodeIndexList *list)
	{
		assert(list);
		base_list = list;
	}
	double solve()
	{
		bestPathLen = std::numeric_limits<double>::max();
		NodeIndexList *usedList;
		usedList = new NodeIndexList(base_list, base_list->getSize());
		usedList->append(base_list->getNode(0)->index);
		solveSub(usedList, 0);
		std::cout << "best: " << bestPathLen << std::endl;
		return bestPathLen;
	}
	double solveSub(NodeIndexList *usedList, double prevPathLen)
	{
		if(prevPathLen > bestPathLen) return 0;
		if(usedList->getSize() >= base_list->getSize()){
			prevPathLen += usedList->getDistanceBetweenFirstAndLast();
			if(bestPathLen > prevPathLen){
				bestPathLen = prevPathLen;
				best_list = new NodeIndexList(usedList, usedList->getSize());
				best_list->appendAll(usedList);
			}
			return 0;
		}
		for(int i = 0; i < base_list->getSize(); i++){
			int index = base_list->getNode(i)->index;
			if(usedList->hasIndex(index)) continue;
			usedList->append(index);
			double newPathLen = prevPathLen + usedList->getLastTwoDistance();
			solveSub(usedList, newPathLen);
			usedList->pop();
		}
		return 0;
	}
	void saveSolution(std::string filename)
	{
		assert(best_list);
		best_list->saveSolution(filename);
	}
	void dumpSolution(std::ofstream *ofs)
	{
		assert(best_list);
		assert(ofs);
		best_list->dumpAllIndex(ofs);
	}
};

class DivideSolver
{
	NodeIndexList *base_list;
	int *clustor_id_list;
	//static const int clustor_size = 16;
	Vector2D *clustor_center;
	static const int num_of_clustors = 2;
	static const int brute_force_capacity = 10;
	DivideSolver **childCluster;
	bool isLeaf;
	int clustorID;
	static int generateCount;
public:
	DivideSolver(NodeIndexList *list)
	{
		assert(list);
		base_list = list;
		clustor_id_list = new int[base_list->getSize()];
		//
		//std::cout << "DivideSolver(size: " << base_list->getSize() << std::endl;
		isLeaf = base_list->getSize() <= brute_force_capacity;
		if(!isLeaf){
			clustering();
		}
		clustorID = DivideSolver::generateCount++;
	}
	void clustering()
	{
		//num_of_clustors = (base_list->getSize() + clustor_size - 1) / clustor_size;
		for(int i = 0; i < base_list->getSize(); i++){
			clustor_id_list[i] = i % num_of_clustors;
		}
		clustor_center = new Vector2D[num_of_clustors];
		do {
			updateClusterCenter();
		} while(updateCluster());
		//dumpClustor();
		genChildClustor();
	}
	void updateClusterCenter()
	{
		for(int k = 0; k < num_of_clustors; k++){
			clustor_center[k] = getClusterCenter(k);
			//std::cout << "cluster center " << k << " ";
			//clustor_center[k].print();
		}
		//std::cout << "-----------" << std::endl;
	}
	int updateCluster()
	{
		int updateCount = 0;
		int nearest_cluster_id;
		double nearest_distance;
		for(int i = 0; i < base_list->getSize(); i++){
			nearest_cluster_id = 0;
			Node *n = base_list->getNode(i);
			nearest_distance = 
				n->calcDistanceTo(&clustor_center[0]);
			for(int k = 1; k < num_of_clustors; k++){
				double distance = n->calcDistanceTo(&clustor_center[k]);
				if(distance < nearest_distance){
					nearest_distance = distance;
					nearest_cluster_id = k;
				}
			}
			if(clustor_id_list[i] != nearest_cluster_id){
				updateCount++;
				clustor_id_list[i] = nearest_cluster_id;
			}
		}
		return updateCount;
	}

	Vector2D getClusterCenter(int cluster_id)
	{
		Vector2D center;
		center.x = 0;
		center.y = 0;
		int count = 0;
		for(int i = 0; i < base_list->getSize(); i++){
			if(clustor_id_list[i] == cluster_id){
				count++;
				center.x += base_list->getNode(i)->x;
				center.y += base_list->getNode(i)->y;
			}
		}
		if(count == 0){
			return clustor_center[cluster_id];
		}
		center.x /= count;
		center.y /= count;
		return center;
	}

	void dumpClustor()
	{
		dumpClustor(0);
	}
	void dumpClustor(int indentDepth){
		if(isLeaf){
			std::cout << "---- Clustor " << clustorID << " ----" << std::endl;
			base_list->printAll();
			return;
		}
		for(int k = 0; k < num_of_clustors; k++){
			childCluster[k]->dumpClustor(indentDepth + 1);
		}	
	}
	void dumpSolution(std::ofstream *ofs)
	{
		if(isLeaf){
			BruteForceSolver bSolver(base_list);
			bSolver.solve();
			bSolver.dumpSolution(ofs);
			return;
		}
		for(int k = 0; k < num_of_clustors; k++){
			childCluster[k]->dumpSolution(ofs);
		}	
	}
	void genChildClustor()
	{
		childCluster = new DivideSolver *[num_of_clustors];
		
		NodeIndexList *nodeIndexList;
		for(int k = 0; k < num_of_clustors; k++){
			nodeIndexList = new NodeIndexList(base_list, base_list->getSize());
			
			for(int i = 0; i < base_list->getSize(); i++){
				if(clustor_id_list[i] == k){
					nodeIndexList->append(this->base_list->getNode(i)->index);
				}
			}
			childCluster[k] = new DivideSolver(nodeIndexList);
		}
	}
};
int DivideSolver::generateCount;
/*
class RoughSolver
{
	public:
	NodeIndexList *base_list;
	double bestPathLen;
	NodeIndexList *best_list;
	static const int updateLimit = 10;
	int updateCount;
	BruteForceSolver(NodeIndexList *list)
	{
		assert(list);
		base_list = list;
	}
	double solve()
	{
		bestPathLen = std::numeric_limits<double>::max();
		NodeIndexList *usedList;
		usedList = new NodeIndexList(base_list, base_list->getSize());
		usedList->append(0);
		updateCount = 0;
		return solveSub(usedList, 0);
	}
	double solveSub(NodeIndexList *usedList, double prevPathLen)
	{
		if(updateCount >= updateLimit) return 0;
		if(prevPathLen > bestPathLen) return 0;
		if(usedList->getSize() >= base_list->getSize()){
			prevPathLen += usedList->getDistanceBetweenFirstAndLast();
			if(bestPathLen > prevPathLen){
				bestPathLen = prevPathLen;
				best_list = new NodeIndexList(usedList, usedList->getSize());
				best_list->appendAll(usedList);
				//
				std::cout << "---- " << prevPathLen << std::endl;
				//best_list->printAll();
				updateCount++;
			}
			return 0;
		}
		for(int i = 0; i < base_list->getSize(); i++){
			int index = base_list->getNode(i)->index;
			if(usedList->hasIndex(i)) continue;
			usedList->append(index);
			double newPathLen = prevPathLen + usedList->getLastTwoDistance();
			solveSub(usedList, newPathLen);
			usedList->pop();
		}
		return 0;
	}
	void saveSolution(std::string filename)
	{
		assert(best_list);
		best_list->saveSolution(filename);
	}
};
*/
/*
class Path
{
private:
	Node **list;
	int size;
	int max_size;
public:
	Path(int max)
	{
		list = new Node *[max];
		max_size = max;
		size = 0;
	}
	int getSize()
	{
		return size;
	}
	bool canAppend()
	{
		return size < max_size;
	}
	void append(Node *node)
	{
		assert(!node->used);
		assert(canAppend());
		list[size] = node;
		size++;
		node->used = true;
	}
	Node *pop()
	{
		assert(size > 0);
		size--;
		Node *node = list[size];
		list[size] = nullptr;
		return node;
	}
	
	Node *findAndAppendNextNode(NodeList *list)
	{
		assert(list);
		Node *bestNode = nullptr;
		for(int i = 0; i < list->getSize(); i++){
			Node *node = &list->list[i];
			assert(node);
			if(node->used) continue;
			node->path_index = -1;
			node->distance_to_path = std::numeric_limits<double>::max();
			for(int k = 0; k < size; k++){
				NodePair pair = getNodePairOfPath(k);
				double distance = pair.calcDistanceFromThisLineSegment(node);
				if(distance < node->distance_to_path){
					node->distance_to_path = distance;
					node->path_index = k;
				}
			}
		}
		for(int i = 0; i < list->getSize(); i++){
			Node *node = &list->list[i];
			if(node->used) continue;
			if(!bestNode || bestNode->distance_to_path > node->distance_to_path){
				bestNode = node;
			}
		}
		if(bestNode){
			appendNodeBetweenPath(bestNode);
		}
		return bestNode;
	}
	void appendNodeBetweenPath(Node *node){
		assert(node);
		assert(size < max_size);
		for(int i = size; i > node->path_index; i--){
			list[i] = list[i - 1];
		}
		list[node->path_index] = node;
		node->used = true;
		size++;
	}
	NodePair getNodePairOfPath(int index)
	{
		assert(0 <= index && index < size);
		return NodePair(list[index], list[(size + index - 1) % size]);
	}
	void saveSolution(std::string filename)
	{
		std::ofstream ofs;
		ofs.open(filename, std::ios::out);
		ofs << "index" << std::endl;
		std::cout << "index" << std::endl;
		for(int i = 0; i < size; i++){
			ofs << list[i]->index << std::endl;
			//std::cout << list[i]->index << std::endl;
		}
	}
};
*/
void TestVector2D()
{
	Vector2D p, q, t;

	p.set(0, 0);
	q.set(0, 0);
	t.set(0, 0);
	assert(t.calcDistanceFromLineSegmentPQ(&p, &q) - 0 < 1e6);

	p.set(0, 0);
	q.set(0, 0);
	t.set(1, 1);
	assert(t.calcDistanceFromLineSegmentPQ(&p, &q) - sqrt(2) < 1e6);

	p.set(0, 0);
	q.set(0, 0);
	t.set(-1, -2);
	assert(t.calcDistanceFromLineSegmentPQ(&p, &q) - sqrt(5) < 1e6);

	p.set(2, 0);
	q.set(0, 2);
	t.set(0, 0);
	assert(t.calcDistanceFromLineSegmentPQ(&p, &q) - sqrt(2) < 1e6);

	p.set(2, 0);
	q.set(0, 2);
	t.set(-2, 4);
	assert(t.calcDistanceFromLineSegmentPQ(&p, &q) - sqrt(2) < 1e6);

	p.set(2, 0);
	q.set(0, 2);
	t.set(4, -2);
	assert(t.calcDistanceFromLineSegmentPQ(&p, &q) - sqrt(2) < 1e6);
	
	std::cout << __FUNCTION__ << ": Test passed." << std::endl;
}

int main(int argc, char *argv[])
{
	TestVector2D();
	//
	std::ifstream ifs;
	if(argc < 3){
		std::cout << "usage: solver <input.csv> <output.csv>" << std::endl;
		return 1;
	}
	//
	NodeList nodeList(MAX_NODE_SIZE);
	nodeList.readFromFile(argv[1]);
	//
	NodeIndexList nodeIndexList(&nodeList, nodeList.getSize());
	nodeIndexList.appendAll(&nodeList);
	//
	DivideSolver dSolver(&nodeIndexList);
	dSolver.dumpClustor();
	std::ofstream ofs;
	ofs.open(argv[2], std::ios::out);
	ofs << "index" << std::endl;
	dSolver.dumpSolution(&ofs);
	//
	return 0;
}
