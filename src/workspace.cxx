#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>

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
  dt= 0.01;
  max_speed = 20.0;
  max_force = 80.0;
  time = 0.,

  this->init();}

Workspace::Workspace(size_t nAgents,
             Real wc, Real wa, Real ws,
             Real rc, Real ra, Real rs) :
             na(nAgents), dt(.05), time(0),
             wCohesion(wc), wAlignment(wa), wSeparation(ws),
             rCohesion(rc), rAlignment(ra), rSeparation(rs),
             max_speed(20.), max_force(80.)
{ this->init(); }

void Workspace::init(){
  lx = 800.0;
  ly = 800.0;
  lz = 800.0;

  padding = 0.02 * lx;
  // Random generator seed
  srand48(std::time(0));

#pragma omp parallel
{
    int k = 0;

    int tid = omp_get_thread_num();
    int max = omp_get_max_threads();

    // Initialize agents
    // This loop may be quite expensive due to random number generation
    for(k = na * (int)(tid/max); k< na * (int)(tid+1/max); k++){
        // Create random position
        //Vector position(lx*(0.02 + drand48()), ly*(0.02 + drand48()), lz*(0.02 + drand48()));
        Vector position(lx*(0.02 + drand48()), ly*(0.02 + drand48()), lz*(0.02 + drand48()));
        Vector velocity(160 * (drand48() - 0.5), 160*(drand48()- 0.5), 160*(drand48() - 0.5));

        // Create random velocity
        agents.push_back(Agent(position, velocity, Zeros()));
        agents.back().max_force = max_force;
        agents.back().max_speed = max_speed;
        agents.back().ra = rAlignment;
        agents.back().rc = rCohesion;
        agents.back().rs = rSeparation;
    }

}

}

void Workspace::move()
{
    // Compute forces applied on specific agent
    for(size_t k = 0; k< na; k++){
      agents[k].compute_force(agents, k, rCohesion);

      agents[k].direction = agents[k].cohesion*wCohesion
        + agents[k].alignment*wAlignment
        + agents[k].separation*wSeparation;
    }

#pragma omp parallel
{
    int k = 0;

    int tid = omp_get_thread_num();
    int max = omp_get_max_threads();

    // Time integraion using euler method
    for(k = na * (int)(tid/max); k< na * (int)(tid+1/max); k++){
      agents[k].velocity += dt*agents[k].direction;

      double speed = agents[k].velocity.norm()/max_speed;
      if (speed > 1. && speed>0.) {
        agents[k].velocity /= speed;
      }
      agents[k].position += dt*agents[k].velocity;

      if(agents[k].position.x <40)
        agents[k].position.x = lx - 40;
      if(agents[k].position.x >lx - 40)
        agents[k].position.x = 40;
      if(agents[k].position.y <40)
        agents[k].position.y = ly - 40;
      if(agents[k].position.y >ly - 40)
        agents[k].position.y = 40;
      if(agents[k].position.z <40)
        agents[k].position.z = lz - 40;
      if(agents[k].position.z >lz - 40)
        agents[k].position.z = 40;
    }

}

}

void Workspace::simulate(int nsteps) {
  // store initial positions
    save(0);

    // perform nsteps time steps of the simulation
    int step = 0;
    while (step++ < nsteps) {
        fprintf(stderr, "steps %i", step);
        this->move();
        // store every 20 steps
        if (step%20 == 0) save(step);
    }
}

void Workspace::save(int stepid) {
  std::ofstream myfile;

  myfile.open("boids.xyz", stepid==0 ? std::ios::out : std::ios::app);

    myfile << std::endl;
    myfile << na << std::endl;
    for (size_t p=0; p<na; p++)
        myfile << "B " << agents[p].position;

    myfile.close();
  }
