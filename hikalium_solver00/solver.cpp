#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <assert.h>

#define MAX_NODE_SIZE	3000

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
	int index, flg;
};

class NodePair
{
public:
	Node *p, *q;
	double calcDistance()
	{
		return p->calcDistanceTo(q);
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
	Node *list;
	int size;
	int max_size;
public:
	NodeList(int max)
	{
		list = new Node[max];
		max_size = max;
		size = 0;
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
};

class Path
{
private:
	int *list;
	int size;
	int max_size;
public:
	Path(int max)
	{
		list = new int[max];
		max_size = max;
		size = 0;
	}
	int append(int nodeIndex)
	{
		if(size >= max_size) return 1;
		list[size] = nodeIndex;
		size++;
		return 0;
	}
	void saveSolution(std::string filename)
	{
		std::ofstream ofs;
		ofs.open(filename, std::ios::out);
		ofs << "index" << std::endl;
		for(int i = 0; i < size; i++){
			ofs << list[i] << std::endl;
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
}

int main(int argc, char *argv[])
{
	TestVector2D();
	std::ifstream ifs;
	if(argc < 3){
		std::cout << "usage: solver <input.csv> <output.csv>" << std::endl;
		return 1;
	}
	ifs.open(argv[1]);
	if(!ifs){
		std::cout << "input file not found." << std::endl;
		return 1;
	}

	NodeList nodeList(MAX_NODE_SIZE);
	Path path(MAX_NODE_SIZE);
	std::string line;
	std::getline(ifs, line);	// skip first line
	while(!ifs.eof()){
		std::string token;
		std::getline(ifs, line);
		std::istringstream stream(line);
		std::string xStr, yStr;
		getline(stream, xStr,',');
		getline(stream, yStr,',');
		try{
			nodeList.append(std::stod(xStr), std::stod(yStr));
		} catch(std::invalid_argument){
			break;
		}
	}
	nodeList.printAll();
	NodePair pair = nodeList.getMostDistantPair();
	pair.print();
	path.append(pair.p->index);
	path.append(pair.q->index);
	path.saveSolution(argv[2]);
	return 0;
}
