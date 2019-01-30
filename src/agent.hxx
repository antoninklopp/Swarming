#ifndef AGENT_HXX
#define AGENT_HXX

#include "types.hxx"
#include "vector.hxx"
#include <vector>
#include "workspace.hxx"

using namespace std;

typedef enum
{
	prey,
	predator,
	active,
	wall
} AgentType;

class Agent
{
  public:

	static int NUMBER_BOIDS;
	static void reinit_number(){Agent::NUMBER_BOIDS=0;}
	Vector position;
	Vector velocity;
	Vector direction;

	Vector cohesion;
	Vector separation;
	Vector alignment;

	double max_speed;
	double max_force;

	// Distance of influence
	double rc, rs, ra;

	bool move;

	int number;

	Agent();

	Agent(const Vector &pos, const Vector &vel, const Vector &dir);

	void compute_force(vector<vector<Agent> > &agent_list, size_t index_list, size_t index_x,
		size_t index_y, size_t index_z, Real lx, Real ly, Real lz, double rad);

	size_t find_closest(vector<Agent> &agent_list, size_t index);
};
#endif
