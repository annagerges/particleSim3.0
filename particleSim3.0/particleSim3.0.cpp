// Executes particle simulation where a user defined amount of particles (100-1000) move around indefinitely. 
// Uses vectors and hashmap based spatial partitioning to store particles and optimize collision detection. Euclidean approximation is used to simulate movement over a fixed time step.

#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <limits>
#include <fstream>
#include <unordered_map>
#include "Particles.h"

using namespace std;

//timestep
const float dt = 0.001f;

int main(int argc, char* argv[])
{
    fstream file("particleInfo.csv",ios::out);

    int nP, row, col, stepCount=0, cellKey;
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
    vector<Particles>particles(nP);

    //sets up random number generator for particle position and velocity
    uniform_real_distribution<float>randPos(7, 799);
    uniform_real_distribution<float>randVelo(1, 30);

    unordered_map <int, vector<Particles*>> hash;

    for (int index = 0; index < nP; index++) {

        //Randomly assigns x and y to be from 1 to 799 because having the user decide would be tedious
        particles[index].setX(randPos(myEngine));
        particles[index].setY(randPos(myEngine));

        //Randomly assigns vx and vy to be from 1 to 30
        particles[index].setVx(randVelo(myEngine));
        particles[index].setVy(randVelo(myEngine));

        //finds cell coordinates based on particle position
        row = particles[index].getY() / width;
        col = particles[index].getX() / width;

        if (row < 0) {
            row = 0;
        }
        else if (row >= nBox) {
            row = nBox - 1;
        }
        if (col < 0) {
            col = 0;
        }
        else if (col >= nBox) {
            col = nBox - 1;
        }

        particles[index].setRow(row);
        particles[index].setCol(col);

        //cell (any paritcle w same row and col will have the same cell key)
        cellKey = (row * nBox) + col;

        hash[cellKey].push_back(&particles[index]);

    }

    
    if (argc > 2) {
        //set k assuming all of the particles are statically laying on the spring and compressing 0.2m. Mulitply by 4 so that the particles are springy.
        s.setK(((nP * 9.8 * particles[0].getMass()) / 0.2) * 4);
        cout << "Too many arguments entered, so k was dynamically allocated to: " << s.getK();
    }
    else if (argc == 2) {
        double k = stof(argv[1]);
        if (k <=0) {
            cout << "INVALID ARGUMENT. K is too small.";
            return 1;
        }
        else {
            //customize k
            s.setK(k);
        }


        cout << "Argument detected! k dynamically allocated to: " << s.getK() << "\n";
    }
    else {
        // If they clicked the normal VS button (argc == 1) or passed too many arguments (argc > 2)
        cout << "Error: You must provide a command-line argument to run this simulation.\n";
        return 1;
    }

    //to run the program locally
    //s.setK(((nP * 9.8 * particles[0].getMass()) / 0.2) * 4);
    
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

            //update position, wall collision checks, and clear and update the hashmap 
            updatePos(particles, s);

            //debugging purposes
             //cout << particles[0].getY() << endl;

            wallCollis(particles);
            hash.clear();
            fix(particles, hash, width, nBox);
            particleCollis(hash);

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
