#include "agent.hxx"
#include <initializer_list>
#include <algorithm>

int Agent::NUMBER_BOIDS;


Agent::Agent()
{
	move = false;
	number = NUMBER_BOIDS;
	NUMBER_BOIDS ++;
}

Agent::Agent(const Vector &pos, const Vector &vel, const Vector &dir) : Agent()
{
	position = pos;
	velocity = vel;
	direction = dir;
	max_speed = 80.0;
	max_force = 20.0;
}

void Agent::compute_force(vector<vector<Agent> > &agent_list, size_t index_list, size_t index_x, size_t index_y, size_t index_z, Real lx, Real ly, Real lz, double rad)
{

	cohesion = Zeros();
	alignment = Zeros();
	separation = Zeros();

	int size_vec_x = (int)(lx / PADDING_GRID);
	int size_vec_y = (int)(ly / PADDING_GRID);
	int size_vec_z = (int)(lz / PADDING_GRID);

	move = false;

	int max_dist = std::max({rs, rc, ra});

	int count_c = 0, count_s = 0, count_a = 0;
	for (int x = max((int)(index_x - (int)(max_dist/PADDING_GRID) - 1), 0); (x < (index_x + (int)(max_dist/PADDING_GRID)) + 1) && (x < size_vec_x); x++){
		for (int y = max((int)(index_y - (int)(max_dist/PADDING_GRID)-1), 0); (y < (index_y + (int)(max_dist/PADDING_GRID)) + 1) && (y<(int)(ly/PADDING_GRID)); y++){
			for (int z = max((int)(index_z - (int)(max_dist/PADDING_GRID)-1), 0); (z < (index_z + (int)(max_dist/PADDING_GRID))+1) && (z<(int)(lz/PADDING_GRID)); z++){
				int k = x * size_vec_y * size_vec_z + y * size_vec_z + z;
				vector<Agent> current_list = agent_list[k];
				for (size_t i = 0; i < current_list.size(); i++)
				{
					Real dist = (this->position - current_list[i].position).norm();
					bool checkNotCurrent = (i==index_list) && (x==index_x) && (y == index_y) && (z == index_z);
					if (!checkNotCurrent && dist < rs && dist > 0.)
					{
						separation += (this->position - current_list[i].position).normalized() / dist;
						++count_s;
					}
					if (!checkNotCurrent && dist < ra)
					{
						alignment += current_list[i].velocity;
						++count_a;
					}
					if (!checkNotCurrent && dist < rc)
					{
						cohesion += current_list[i].position;
						++count_c;
					}
				}
			}
		}
	}

	// if (number == 2000){
	// 	cout << "influence " << count_s << " " << count_a << " " << count_c << endl;
	// 	cout << "velocity" << velocity.x << " " << velocity.y << endl;
	// }

	// Compute separation contribution
	if (count_s > 0)
	{
		separation.normalize();
		separation *= max_force;
		separation -= velocity;
		if (separation.norm() > max_force)
		{
			separation.normalize();
			separation *= max_force;
		}
	}
	// Compute alignment contribution
	if (count_a > 0)
	{
		if (alignment.norm() > 0.0)
		{
			alignment.normalize();
			alignment *= max_force;
		}
		alignment -= velocity;
		if (alignment.norm() > max_force)
		{
			alignment.normalize();
			alignment *= max_force;
		}
	}
	// Compute cohesion contribution
	if (count_c > 0)
	{
		// Compute center of gravity
		cohesion /= count_c;

		// Direction of displacement
		cohesion -= position;

		cohesion.normalize();
		cohesion *= max_speed;
		cohesion -= velocity;
		if (cohesion.norm() > max_force)
		{
			cohesion.normalize();
			cohesion *= max_force;
		}
	}
}

size_t Agent::find_closest(vector<Agent> &agent_list, size_t index)
{
	size_t closest_agent = index;
	double min_dist = 1000;

	double dist;

	for (size_t i = 0; i < agent_list.size(); i++)
	{
		if (i != index)
		{
			dist = (this->position - agent_list[i].position).norm();

			if (dist < min_dist)
			{
				min_dist = dist;
				closest_agent = i;
			}
		}
	}
	cout << "closest";
	return closest_agent;
}
