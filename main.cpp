using namespace std;
#include <vector>
#include <stdio.h>
#include <string.h>
#include <ctime>
#include "engine.h"

//use seed from time to generate random values
std::mt19937 rng(time(0));

int main(int argc, char** argv)
{
    string line;
    vector<Particle> particles;
    int N, L, r, S, i = 0;
    double x, y, vx, vy;
    char mode_buf[6];

    scanf("%d", &N);
    scanf("%d", &L);
    scanf("%d", &r);
    scanf("%d", &S);
    scanf("%5s", mode_buf);

    bool printing = (strcmp(mode_buf, "print") == 0);

    //Make N particles
    for (i = 0; i < N; i++) {
        Particle newpart(-1, 0, 0, 0, 0);
        particles.push_back(newpart);
    }

    i = 0; //added
    while (scanf("%d %lf %lf %lf %lf", &i, &x, &y, &vx, &vy) != EOF) {
        particles[i].setN(i);
        particles[i].setX(x);
        particles[i].setY(y);
        particles[i].setVx(vx);
        particles[i].setVy(vy);
    }

    //if values were not given in input file
    if(particles[0].getN() == -1){
        //generate a distribution for random values in given range
        std::uniform_int_distribution<int> dist_pos(r, L-r);

        std::uniform_int_distribution<int> dist_v(L/(8*r), L/4);

        //Generate particles and store into data structure
        for(int j = 0; j < N; j++){

            x=dist_pos(rng);
            y=dist_pos(rng);
            vx=dist_v(rng);
            vy=dist_v(rng);

            particles[j].setN(j);
            particles[j].setX(x);
            particles[j].setY(y);
            particles[j].setVx(vx);
            particles[j].setVy(vy);
        }
    }

    //start simulation
    Engine eng(particles, N, L, r, S, printing);

    return 0;
}

//srand(0);
