#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

#include "agent.hxx"
#include "vector.hxx"
#include "workspace.hxx"
#include <omp.h>

Workspace::Workspace(ArgumentParser &parser)
{

	na = parser("agents").asInt();

	wCohesion = parser("wc").asDouble();
	wAlignment = parser("wa").asDouble();
	wSeparation = parser("ws").asDouble();

	rCohesion = parser("rc").asDouble();
	rAlignment = parser("ra").asDouble();
	rSeparation = parser("rs").asDouble();
	dt = 0.01;
	max_speed = 20.0;
	max_force = 80.0;
	cout << "size" << agents.size() << endl;
	time = 0.,

	this->init();
}

Workspace::Workspace(size_t nAgents,
					 Real wc, Real wa, Real ws,
					 Real rc, Real ra, Real rs) : na(nAgents), dt(.05), time(0),
												  wCohesion(wc), wAlignment(wa), wSeparation(ws),
												  rCohesion(rc), rAlignment(ra), rSeparation(rs),
												  max_speed(20.), max_force(80.)
{
	this->init();
}

void Workspace::init()
{
	lx = 800.0;
	ly = 800.0;
	lz = 800.0;

	padding = 0.02 * lx;
	// Random generator seed
	srand48(std::time(0));

	int k = 0;
	int size_vec_x = (int)(lx / PADDING_GRID);
	int size_vec_y = (int)(ly / PADDING_GRID);
	int size_vec_z = (int)(lz / PADDING_GRID);
	
	agents.resize(size_vec_x * size_vec_y * size_vec_z);

#pragma omp parallel private(k) shared(size_vec_z, size_vec_y, size_vec_x)
	{

		int tid = omp_get_thread_num();
		int max = omp_get_max_threads();

		cout << "threads " << tid << endl;

		// Initialize agents
		// This loop may be quite expensive due to random number generation
		for (int div_x = 0; div_x < PADDING_CORE; div_x ++){
			for (int div_y = 0; div_y < PADDING_CORE; div_y ++){
				for (int div_z = 0; div_z < PADDING_CORE; div_z ++){
					if (div_x * PADDING_CORE * PADDING_CORE + div_y * PADDING_CORE + div_z % max != tid){
						continue;
					}
					for (size_t x = (int)((div_x * size_vec_x)/PADDING_CORE); x < (int)(((div_x+1) * size_vec_x)/PADDING_CORE); x++){
						for (size_t y = (int)((div_y * size_vec_y)/PADDING_CORE); y < (int)(((div_y+1) * size_vec_y)/PADDING_CORE); y++){
							for (size_t z = (int)((div_z * size_vec_z)/PADDING_CORE); z < (int)(((div_z+1) * size_vec_z)/PADDING_CORE); z++){
								k = x * PADDING_GRID * PADDING_GRID + y * PADDING_GRID + z;
								for (size_t i = 0; i < agents[k].size(); i++){
									cout << "k" << k  << "size_vec" << agents.size() << endl;
									// Create random position
									//Vector position(lx*(0.02 + drand48()), ly*(0.02 + drand48()), lz*(0.02 + drand48()));
									Vector position(lx * (0.02 + drand48()), ly * (0.02 + drand48()), lz * (0.02 + drand48()));
									Vector velocity(160 * (drand48() - 0.5), 160 * (drand48() - 0.5), 160 * (drand48() - 0.5));

									// Create random velocity
									agents[k].push_back(Agent(position, velocity, Zeros()));
									agents[k].back().max_force = max_force;
									agents[k].back().max_speed = max_speed;
									agents[k].back().ra = rAlignment;
									agents[k].back().rc = rCohesion;
									agents[k].back().rs = rSeparation;
								}
							}
						}
					}
				}
			}
		}
	}
}

void Workspace::move()
{

	int k = 0;
	int size_vec_x = (int)(lx / PADDING_GRID);
	int size_vec_y = (int)(ly / PADDING_GRID);
	int size_vec_z = (int)(lz / PADDING_GRID);

#pragma omp parallel private(k) shared(size_vec_z, size_vec_y, size_vec_x)
	{

		int tid = omp_get_thread_num();
		int max = omp_get_max_threads();

		for (int div_x = 0; div_x < PADDING_CORE; div_x ++){
			for (int div_y = 0; div_y < PADDING_CORE; div_y ++){
				for (int div_z = 0; div_z < PADDING_CORE; div_z ++){
					if (div_x * PADDING_CORE * PADDING_CORE + div_y * PADDING_CORE + div_z % max != tid){
						continue;
					}
					for (size_t x = (int)((div_x * size_vec_x)/PADDING_CORE); x < (int)(((div_x+1) * size_vec_x)/PADDING_CORE); x++){
						for (size_t y = (int)((div_y * size_vec_y)/PADDING_CORE); y < (int)(((div_y+1) * size_vec_y)/PADDING_CORE); y++){
							for (size_t z = (int)((div_z * size_vec_z)/PADDING_CORE); z < (int)(((div_z+1) * size_vec_z)/PADDING_CORE); z++){
								k = x * PADDING_GRID * PADDING_GRID + y * PADDING_GRID + z;
								for (size_t i = 0; i < agents[k].size(); i++){
									agents[k][i].compute_force(agents, i, x, y, z, lx, ly, lz, rCohesion);

									#pragma omp critical
									agents[k][i].direction = agents[k][i].cohesion * wCohesion +
										agents[k][i].alignment * wAlignment + agents[k][i].separation * wSeparation;
								}
							}
						}
					}
				}
			}
		}

		// Time integraion using euler method
		// for (k = (int)(na * tid / max); k < (int)(na * (tid + 1) / max); k++)
		// {
		// 	agents[k].compute_force(agents, k, rCohesion);
		//
		// 	#pragma omp critical
		// 	agents[k].direction = agents[k].cohesion * wCohesion + agents[k].alignment * wAlignment + agents[k].separation * wSeparation;
		// }
	}

#pragma omp parallel private(k)
	{

		int tid = omp_get_thread_num();
		int max = omp_get_max_threads();

		for (k = (int)(tid * PADDING_GRID * PADDING_GRID * PADDING_GRID/max);
			k <(int)((tid+1) * PADDING_GRID * PADDING_GRID * PADDING_GRID/max); k++){

			// Time integraion using euler method
			for (int i = 0; i < agents[k].size(); i++)
			{
				agents[k][i].velocity += dt * agents[k][i].direction;

				double speed = agents[k][i].velocity.norm() / max_speed;
				if (speed > 1. && speed > 0.)
				{
					agents[k][i].velocity /= speed;
				}
				agents[k][i].position += dt * agents[k][i].velocity;

				if (agents[k][i].position.x < 40)
					agents[k][i].position.x = lx - 40;
				if (agents[k][i].position.x > lx - 40)
					agents[k][i].position.x = 40;
				if (agents[k][i].position.y < 40)
					agents[k][i].position.y = ly - 40;
				if (agents[k][i].position.y > ly - 40)
					agents[k][i].position.y = 40;
				if (agents[k][i].position.z < 40)
					agents[k][i].position.z = lz - 40;
				if (agents[k][i].position.z > lz - 40)
					agents[k][i].position.z = 40;
			}
		}
	}
}

void Workspace::simulate(int nsteps)
{
	// store initial positions
	save(0);

	// perform nsteps time steps of the simulation
	int step = 0;
	while (step++ < nsteps)
	{
		this->move();
		cout << "steps " << step << endl;
		// store every 20 steps
		// if (step % 20 == 0)
		// 	save(step);
	}
	cout << "FINISHED" << endl;
}

void Workspace::save(int stepid)
{
	std::ofstream myfile;

	myfile.open("boids.xyz", stepid == 0 ? std::ios::out : std::ios::app);

	myfile << std::endl;
	myfile << na << std::endl;
	for (size_t p = 0; p < agents.size(); p++)
		for (size_t i = 0; i < agents[p].size(); i++)
			myfile << "B " << agents[p][i].position;

	myfile.close();
}
