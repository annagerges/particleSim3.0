// Executes particle simulation where a user defined amount of particles (100-1000) move around indefinitely. 
// Uses vectors and grid-based spatial partitioning to store particles and optimize collision detection. Euclidean approximation is used to simulate movement over a fixed time step.

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <limits>
#include <fstream>
#include "Particles.h"

using namespace std;

//timestep and width of each grid cell
const float dt = 0.001f;

int main()
{
    fstream file("particleInfo.csv",ios::out);

    int nP, row, col, stepCount=0;
    float accumulator, totalTime=0;
    Spring s;


    random_device myEngine;

    //prompts user to enter valid amount of particles
    cout << "How many particles do you want(100-1000): ";
    cin >> nP;

    int nBox = max(5, (int)sqrt(nP / 5));
    int width = 800 / nBox;

    while (nP < 100 || nP>1000) {
        cout << "\nInvalid Input. Enter a valid amount: ";
        cin >> nP;
    }


    //creates a vector of particles with the valid size the user specified to keep track of every particle
    vector<Particles>particles;
    particles.reserve(nP);

    //the grid (3d vector of Particle pointers) to optimize the program
    vector<vector<vector<Particles*>>>grid(nBox, vector<vector<Particles*>>(nBox));

    // Pre-reserve space in each grid cell to avoid reallocations
    int estimatedPerCell = nP / (nBox * nBox) + 2;  // +2 for buffer
    for (int i = 0; i < nBox; i++) {
        for (int j = 0; j < nBox; j++) {
            grid[i][j].reserve(estimatedPerCell);
        }
    }

    if (argc > 1) {
        //set k assuming all of the particles are statically laying on the spring and compressing 0.2m. Mulitply by 4 so that the particles are springy.
        s.setK(((nP * 9.8 * particles[0].getMass()) / 0.2) * 4);
        cout << "No argument was entered for k, so it was dynamically allocated to: " << s.getK();
    }

    //sets up random number generator for particle position and velocity
    uniform_real_distribution<float>randPos(7, 799);
    uniform_real_distribution<float>randVelo(1, 30);

    for (int index = 0; index < nP; index++) {
        particles.emplace_back();

        //Randomly assigns x and y to be from 1 to 799 because having the user decide would be tedious
        particles[index].setX(randPos(myEngine));
        particles[index].setY(randPos(myEngine));

        //Randomly assigns vx and vy to be from 1 to 30
        particles[index].setVx(randVelo(myEngine));
        particles[index].setVy(randVelo(myEngine));

        //finds grid cell based on particle position
        row = particles[index].getY() / width;
        col = particles[index].getX() / width;

        int nR = grid.size(), nC = grid[0].size();

        if (row < 0) {
            row = 0;
        }
        else if (row >= nR) {
            row = nR - 1;
        }
        if (col < 0) {
            col = 0;
        }
        else if (col >= nC) {
            col = nC - 1;
        }

        particles[index].setRow(row);
        particles[index].setCol(col);

        //puts the particle in the cell
        grid[row][col].push_back(&particles[index]);
    }

    //set k assuming all of the particles are statically laying on the spring and compressing 0.2m. Mulitply by 4 so that the particles are springy.
    s.setK(((nP * 9.8 * particles[0].getMass()) / 0.2) * 4);

    //writing k and num of particles into the file for it to be analyzed using python but not seen
    file << "# nP: " << nP << "\n";
    file << "# k: " << s.getK() << "\n";
    file << "# h: " << s.getHeight() << "\n";

    //file rows
    file << "Particle num,x,y,vx,vy,ay,cellRow,cellCol,time(s)" << "\n";

    //sets accumulator to 0
    accumulator = 0;

    //start frame timer
    auto previousTime = chrono::high_resolution_clock::now();


    //infinite loop
    while (true) {
        //end of frame timer
        auto currentTime = chrono::high_resolution_clock::now();

        //duration of frame in seconds
        chrono::duration<float> frameDur = currentTime - previousTime;

        //assign the start of the next frame to be the time the previous one ended
        previousTime = currentTime;

        //increment accumulator by frame duration
        accumulator += frameDur.count();

        while (accumulator >= dt) {
            totalTime += dt;

            //update position, wall collision checks, and clear and update the grid
            updatePos(particles, s);

            //debugging purposes
             //cout << particles[0].getY() << endl;

            wallCollis(particles);
            clearAndFix(particles, grid, width);

            for (int row = 0; row < grid.size(); row++) {
                for (int col = 0; col < grid[0].size(); col++) {
                    //if a grid cell has more then 1 particle (is active) then check if they are colliding
                    if (grid[row][col].size() > 1) {
                        particleCollis(grid[row][col]);
                    }
                }
            }
            stepCount++;

            //every 10 frames (0.1 seconds) every particle state is logged into csv file for further analysis
            if (stepCount % 10 == 0) {
                csvDump(particles, file,totalTime);
            }

            //decrement accumulator by the timestep
            accumulator -= dt;
        }
    }

    file.close();

    return 0;
}
