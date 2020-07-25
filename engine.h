#ifndef ENGINE_H
#define ENGINE_H
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <map>
#include <vector>
#include <forward_list>

#include "vector.h"
#include "particle.h"

using namespace std;

//used to set time to collision as a high value to indicate no collision
int const NO_VALUE = std::numeric_limits<int>::max();

//type of collision
enum CollType { NOCOLL, LEFTWALL, RIGHTWALL, TOPWALL, BOTTOMWALL, CORNER, PARTCOLL };

//used to return what kind of collision happened with default values if
//a collision did not happen
struct Collision {
	//colliding particle #1 (only one colliding if hitting a wall)
	int partId1 = NO_VALUE;
	//colliding particle #2
	int partId2 = NO_VALUE;
	//tells if we collided with a particle, with which wall
	//or didn't collide
	CollType type = NOCOLL;
	//time it took to collide
	double t = NO_VALUE;
	//wall that particles stop at after colliding
	CollType next1 = NOCOLL;
	CollType next2 = NOCOLL;
};


class Engine
{
public:
	//constructor and deconstructor
	Engine(vector<Particle> particles,
		int N, int L, int r, int S, bool printing);
	~Engine();

	//Used to print info on paticles (velocity, coordinate etc.)
	void print();

	//advancing in the simulation
	void advanceStep();

	static bool earlier(Collision c1, Collision c2);

	//finding all collisions within the given time t
	//info of collision is writen to the provided list
	void allCollisions(forward_list<Collision> &collisions, double t);

private:
	//contains all particles
	//map<int, Particle> particles_;
	vector<Particle> particles_;
	//is filled with collisions
	forward_list<Collision> collisions;

	//amount of particles
	int N_;
	//length and height of square box
	int L_;
	//rdius of particles
	int r_;
	//total amount of steps in simulation
	int S_;
	//"print" to print info on each step, "printf" to print only initial and final info
	bool printing_;

	//moves particles by time t, unless particles are stopped at a wall
	//due to multple collisions
	void advanceParticles(double t);

	//finds all of a particle's particle collisions within t, from their current positions
	//discount duplicates?
	void findPartCols(double t, Particle& part1, forward_list<Collision> &collisions);

	//finds a particle's collision with a wall within time t, from its current position
	Collision getWallCol(double t, Particle& part1);

	//current step
	int curS_ = 0;
};

//finds time it takes for given particles to collide with each other from their current positions
//if it happens within t, collision info is returned
Collision particleCollisionTime(Particle& part1, Particle& part2, double t, int r);

//collides given particles (changes their velocities to what they will be after collision)
void collideParticles(Particle& part1, Particle& part2);

//collides particle with the given wall
void collideWall(Particle& part1, CollType col);


#endif // ENGINE_H
