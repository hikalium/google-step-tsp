#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>

#define MAX_NODE_SIZE	3000

typedef struct {
	double x, y;
	int flg;
} Node;

Node node_list[MAX_NODE_SIZE];
int num_of_nodes = 0;
int tmp_route[MAX_NODE_SIZE];
int min_route[MAX_NODE_SIZE];
double min_cost;
/*
int solveByBruteForce()
{
	int min_cost = 
}
*/
double calcDistance(Node *a, Node *b)
{
	double dx = (a->x - b->x);
	double dy = (a->y - b->y);
	return sqrt(dx * dx + dy * dy);
}

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

void saveSolution(std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename, std::ios::out);
	ofs << "index" << std::endl;
	for(int i = 0; i < num_of_nodes; i++){
		ofs << min_route[i] << std::endl;
	}
}

int main(int argc, char *argv[])
{
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

	int lastFrom = 0;
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
			node_list[num_of_nodes].x = std::stod(xStr);
			node_list[num_of_nodes].y = std::stod(yStr);
			num_of_nodes++;
		} catch(std::invalid_argument){
			break;
		}
	}
	for(int i = 0; i < num_of_nodes; i++){
		std::cout << node_list[i].x << ":" << 
			node_list[i].y << std::endl;
	}
	solveRough();
	saveSolution(argv[2]);
	return 0;
}
