/*
*/
#ifndef WORKSPACE
#define WORKSPACE

#include "parser.hxx"
#include "types.hxx"

const int PADDING_CORE=4;
const float PADDING_GRID=50;

using namespace std;

class Workspace
{
  protected:
	vector<vector<Agent> > agents;
	unsigned int na;

	Real dt;
	int time;

	Real wCohesion, wAlignment, wSeparation;
	Real rCohesion, rAlignment, rSeparation;
	Real maxU;

	Real max_speed;
	Real max_force;

	Real tUpload, tDownload, tCohesion, tAlignment, tSeparation;

	// Size of the domain
	Real lx, ly, lz;

	// Lower bound of the domain
	Real xmin, ymin, zmin;

	// Padding around the domain
	Real padding;

	Real domainsize;
	void init();

  public:
	Workspace(ArgumentParser &parser);

	Workspace(size_t nAgents,
			  Real wc, Real wa, Real ws,
			  Real rc, Real ra, Real rs);

	void move();
	void simulate(int nsteps);
	void save(int stepid);
};

#endif
