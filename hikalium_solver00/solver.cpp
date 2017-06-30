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
public:
	int index, path_index;
	double distance_to_path;
	bool used;
	Node(){
		used = false;
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
	/*
	NodeList getNodeListWithoutNode(NodeList *nodeList, Node *node)
	{
		NodeList *newList = new NodeList(nodeList->getSize() - 1);
		for(int i = 0; i < nodeList->getSize(); i++){
			if(nodeList->list[i] != node){
				newList.append(nodeList->list[i]);
			}
		}
		return NodeList;
	}
	*/
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
	int max_size;
	int size;
public:
	NodeIndexList(NodeList *base_list, int max_size)
	{
		assert(base_list);
		this->base_list = base_list;
		list = new int[max_size];
		this->max_size = max_size;
		this->size = 0;
	}
	int append(int index)
	{
		if(size >= max_size) return 1;
		list[size] = index;
		size++;
		return 0;
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
	void printAll()
	{
		for(int i = 0; i < size; i++){
			base_list->getNode(list[i])->print();
		}
	}
};


/*
class BruteForceSolver
{
	public:
	Path *path;
	double bestPathLen;
	int *bestPath;
	BruteForceSolver(NodeList *list)
	{
		assert(list);
		path = new Path(list->getSize());
		bestPath = new int[list->getSize()];
	}
	double solve()
	{
		bestPathLen = std::numeric_limits<double>::max();
		solveSub(0);
	}
	double solveSub(double prevPathLen)
	{
		if(!path.canAppend()){
			if(prevPathLen < bestPathLen){
				bestPathLen = prevPathLen;
				for(int i = 0; i < list->getSize(); i++){
					
				}
			}
			return;
		}
	}
}
*/
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

/*
int solveByBruteForce()
{
	int min_cost = 
}
*/
/*
void solveRough()
{
	min_cost = 0;
	min_route[0] = 0;
	for(int i = 1; i < num_of_nodes; i++){
		min_route[i] = i;
		min_cost += calcDistance(&node_list[i - 1], &node_list[i]);
	}
	min_cost += calcDistance(&node_list[0], &node_list[num_of_nodes - 1]);
	std::cout << "solveRough: cost = " << min_cost << std::endl;
}
*/
/*
void saveSolution(std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename, std::ios::out);
	ofs << "index" << std::endl;
	for(int i = 0; i < num_of_nodes; i++){
		ofs << min_route[i] << std::endl;
	}
}
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
	std::ifstream ifs;
	if(argc < 3){
		std::cout << "usage: solver <input.csv> <output.csv>" << std::endl;
		return 1;
	}
	//
	NodeList nodeList(MAX_NODE_SIZE);
	nodeList.readFromFile(argv[1]);
	Path path(MAX_NODE_SIZE);
	NodePair pair = nodeList.getMostDistantPair();
	pair.print();
	path.append(pair.p);
	path.append(pair.q);
	while(path.getSize() < nodeList.getSize()){
		path.findAndAppendNextNode(&nodeList);
	}
	path.saveSolution(argv[2]);
	return 0;
}
