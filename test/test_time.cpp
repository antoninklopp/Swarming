
// Agent (particle model)
#include "types.hxx"
#include "parser.hxx"

#include "workspace.cxx"
#include <time.h>

// Main class for running the parallel flocking sim
void test(int number_agents, int steps, float padding)
{
	// Create parser
	ArgumentParser parser;

	// Add options to parser
	parser.addOption("agents", number_agents);
	parser.addOption("steps", steps);
	parser.addOption("wc", 1.0); //7.0);
	parser.addOption("wa", 1.0); // 12.0);
	parser.addOption("ws", 1.5); // 55.0);

	parser.addOption("rc", 90);
	parser.addOption("ra", 90);
	parser.addOption("rs", 25);

	// Create workspace
	Workspace workspace(parser);

	// Launch simulation
	int nSteps = parser("steps").asInt();
	workspace.simulate(nSteps);
}

int main(){
		clock_t tStart = clock();
		struct timespec start, finish;
		double elapsed;

		// for (int steps = 10; steps < 300; steps+=10){
		//     tStart = clock();
		//     test(10000, steps);
		//     printf("Time taken for %i steps: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC, steps);
		// }

		for (int agents = 1000; agents < 30000; agents+=1000){
				elapsed = 0;
				clock_gettime(CLOCK_MONOTONIC, &start);
				test(agents, 20, 50);
				clock_gettime(CLOCK_MONOTONIC, &finish);

				elapsed = (finish.tv_sec - start.tv_sec);
				elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
				printf("Time taken for %i agents: %.5f s\n", elapsed, agents);
		}
}
