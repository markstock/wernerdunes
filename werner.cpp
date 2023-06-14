//
// werner - a parallel Werner algorithm for modeling sand dunes
//
// (c)2023 Mark J. Stock <markjstock@gmail.com>
//

#include "inout.h"
#include "CLI11.hpp"

#include <Eigen/Core>

#include <cassert>
#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <random>


// begin execution here

int main(int argc, char const *argv[]) {

  std::cout << "werner v0.1\n";

  // process command line args
  CLI::App app{"Model sand dune formation and evolution using Werner's algorithm"};

  // load a dem from a png file - check command line for file name
  std::string demfile = "uniform";
  app.add_option("-g,--ground", demfile, "png DEM for no-sand, unerodable surface");
  std::string sandfile = "uniform";
  app.add_option("-s,--sand", sandfile, "png DEM with sand height");

  // the vertical scale of black-to-white
  int32_t vscale = 4;
  app.add_option("-v,--vert", vscale, "number of vertical steps per horizontal pixel");

  // if demfile is "uniform" or "slope" or "cone" or "random", set the nx, ny here
  size_t nx = 1000;
  app.add_option("-x,--nx", nx, "number of pixels in horizontal direction (if no dem png is given)");
  size_t ny = 1000;
  app.add_option("-y,--ny", ny, "number of pixels in vertical direction (if no dem png is given)");

  // runtime parameters
  int32_t ld = 5;
  app.add_option("-l,--ldist", ld, "distance to skip every step, >=1");
  float pns = 0.4;
  app.add_option("-n,--pns", pns, "probability that a hit on ground sticks, 0..1");
  float ps = 0.6;
  app.add_option("-p,--ps", ps, "probability that a hit on sand sticks, 0..1");
  int32_t aor = 4;
  app.add_option("-r,--repose", ld, "angle of repose in steps, >=1");
  int32_t ashade = 2;
  app.add_option("-d,--shadow", ld, "shadow angle in steps, >=1");

  // random seed, if not set
  uint32_t rseed = 12345;
  app.add_option("--seed", rseed, "random seed, program will use noise generator if this is not set");

  // write out certain arrays
  std::string outsand;
  app.add_option("-o,--out", outsand, "write sand height to png");

  // finally parse
  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }


  //
  // Preliminary work - load in the arrays
  //

  // load the data

  // array of the ground elevations (DTM)
  Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic> ground;

  if (demfile.empty() or demfile == "uniform") {
    std::cout << "Setting elevations to constant value (flat)\n";

    ground.resize(nx,ny);
    ground.setZero();

  } else if (demfile == "slope") {
    std::cout << "Setting elevations to uniform slope\n";

    ground.resize(nx,ny);
    for (size_t i=0; i<nx; ++i) for (size_t j=0; j<ny; ++j) ground(i,j) = i;

  } else if (demfile == "cone") {
    std::cout << "Setting elevations to smooth cone\n";

    ground.resize(nx,ny);
    const float ccdist = std::sqrt(0.5);
    for (size_t i=0; i<nx; ++i) for (size_t j=0; j<ny; ++j) {
      ground(i,j) = vscale * (float)nx * (1.0 - std::sqrt(std::pow(0.5-i/(float)(nx-1), 2) + std::pow(0.5-j/(float)(ny-1), 2))/ccdist);
      if (i%10==0 and j%10==0) std::cout << i << " " << j << " is " << ground(i,j) << "\n";
    }

  } else if (demfile == "random") {
    std::cout << "Setting elevations to random noise\n";

    // since Eigen uses rand() internally, we can set the seed here
    if (app.count("--seed") > 0) {
      // a random seed was given on the command line
      srand((unsigned int) rseed);
    } else {
      //seed based on a low-resolution timer
      srand((unsigned int) time(0));
    }
    ground.resize(nx,ny);
    ground.setRandom();	// sets to [-1..1]
    ground = ground.array().pow(2);		// change to all positive

  } else {
    // read a png to get the elevation
    // check the resolution first
    int hgt, wdt;
    (void) read_png_res (demfile.c_str(), &hgt, &wdt);
    if (wdt > 0) nx = wdt;
    if (hgt > 0) ny = hgt;

    // allocate the space
    ground.resize(nx,ny);
    float** data = allocate_2d_array_f((int)nx, (int)ny);

    // read the first channel into the elevation array, scaled as 0..vscale
    (void) read_png (demfile.c_str(), (int)nx, (int)ny, 0, 0, 0.0, 0,
                     data, 0.0, (float)nx, nullptr, 0.0, 1.0, nullptr, 0.0, 1.0);

    for (size_t i=0; i<nx; ++i) for (size_t j=0; j<ny; ++j) ground(i,j) = data[i][j];

    free_2d_array_f(data);
  }

  // write out corners of matrix
  //std::cout << "Top left corner of elevation matrix:" << std::endl;
  //std::cout << ground.block(0,0,6,6) << std::endl;
  //std::cout << std::endl;
  //std::cout << "Bottom right corner of elevation matrix:" << std::endl;
  //std::cout << ground.block(nx-6,ny-6,6,6) << std::endl;


  // another matrix which would increase the cost of traversal
  // this could be used for rivers, woods, etc.
  Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic> sand;

  if (sandfile.empty() or sandfile == "uniform") {
    std::cout << "Setting sand to constant value (flat)\n";

    sand.resize(nx,ny);
    sand.setZero();
    for (size_t i=0; i<nx; ++i) for (size_t j=0; j<ny; ++j) sand(i,j) = 10*vscale;

  } else {
    // check the resolution first
    int hgt, wdt;
    (void) read_png_res (sandfile.c_str(), &hgt, &wdt);
    if (wdt != nx) exit(1);
    if (hgt != ny) exit(1);

    // allocate the space
    sand.resize(nx,ny);
    float** data = allocate_2d_array_f((int)nx, (int)ny);

    // read the first channel into a temp array, scaled as 0..vscale
    (void) read_png (sandfile.c_str(), (int)nx, (int)ny, 0, 0, 0.0, 0,
                     data, 0.0, (float)nx, nullptr, 0.0, 1.0, nullptr, 0.0, 1.0);

    for (size_t i=0; i<nx; ++i) for (size_t j=0; j<ny; ++j) sand(i,j) = data[i][j];

    free_2d_array_f(data);
  }


  // array to store shadow height in
  Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic> shadow;
  shadow.resize(nx,ny);
  shadow.setZero(nx,ny);

  // calculate the shadow map
  // first set it to sand, as it can never be lower
  shadow = sand;

  // now march over the input (left) edge, as flow moves to the right
  for (int32_t j=0; j<ny; ++j) {
    const int32_t ty = j % ny;
    int32_t shadow_height = sand(0,ty);
    // go across the field twice because we are periodic, and shadows may cross the boundary
    for (int32_t i=1; i<2*nx; ++i) {
      const int32_t tx = i % nx;
      shadow_height -= ashade;
      shadow_height = std::max(shadow_height, sand(tx,ty));
      shadow(tx,ty) = shadow_height;
    }
  }


  // use random device
  std::random_device dev;
  std::mt19937 rng(dev());
  // or use seed, if given
  if (app.count("--seed") > 0) {
    rng.seed(rseed);
  }
  // use these to select a random point in the array
  std::uniform_int_distribution<std::mt19937::result_type> xrand(0,nx-1);
  std::uniform_int_distribution<std::mt19937::result_type> yrand(0,ny-1);
  // and use this for probability checks
  std::uniform_real_distribution<> prob(0.0, 1.0);


  //
  // loop over all given start points
  //

  const size_t n_per_step = 10*nx*ny;
  const size_t max_steps = 1000;

  for (size_t step=0; step<max_steps; ++step) {
    std::cout << "\ntime step " << step << "\n";

    // move many sand grains
    for (size_t ipt=0; ipt<n_per_step; ++ipt) {

      const bool debug = (ipt % 1000000 == 0);
      //const bool debug = true;
      if (debug) std::cout << "  saltating grain " << ipt << "\n";

      // find a seed point
      int32_t px,py;
      bool no_sand = true;
      while (no_sand) {
        px = xrand(rng);
        py = yrand(rng);
        //std::cout << "    checking " << px << " " << py << " for sand\n";
        if (sand(px,py) > ground(px,py)) no_sand = false;
      }
      if (debug) std::cout << "    moving sand from " << px << " " << py << "\n";

      // pick it up
      sand(px,py) -= 1;

      // lower the shadow map if necessary
      if (sand(px,py)+1 == shadow(px,py)) {
        // yes, we were the cell casting the shadow (or were exactly on it)
        if (debug) std::cout << "    lowering shadow\n";
        int32_t shadow_height = sand(px,py);
        while (shadow_height >= sand(px,py)) {
          shadow(px,py) = shadow_height;
          shadow_height -= ashade;
          px = (px + 1) % nx;
        }
        //if (debug and shadow_length > 0) std::cout << "    casts shadow " << shadow_length << "\n";
      }

      // move it
      bool keep_moving = true;
      while (keep_moving) {
        // eventually have random variation of the hop distance
        px = (px + ld) % nx;
        if (debug) std::cout << "    at " << px << " " << py << " ground=" << ground(px,py) << " sand=" << sand(px,py) << " shadow=" << shadow(px,py) << "\n";

        // main sticking logic is here!
        if (shadow(px,py) > sand(px,py)) {
          // we're in the shadow region, always stick
          keep_moving = false;
        } else if (sand(px,py) > ground(px,py)) {
          // there is sand here, do we stick?
          if (prob(rng) < ps) keep_moving = false;
        } else {
          // there is no sand here, do we stick?
          if (prob(rng) < pns) keep_moving = false;
        }
      }
      if (debug) std::cout << "    sticks at " << px << " " << py << "\n";

      // now we have a test position, check vs. angle of repose
      int32_t lowx = px;
      int32_t lowy = py;
      bool keep_slumping = true;
      // repeat the slump operation until it can't move any more
      while (keep_slumping) {
        int32_t lowest_slope = 0;
        for (int32_t i=px-1; i<px+2; ++i) {
          const int32_t tx = i % nx;
          for (int32_t j=py-1; j<py+2; ++j) {
            const int32_t ty = j % ny;
            const int32_t slope = sand(tx,ty) - sand(px,py);
            if (slope < lowest_slope) {
              lowx = tx;
              lowy = ty;
              lowest_slope = slope;
            }
          }
        }
        if (debug) std::cout << "    lowest slope is " << lowest_slope << "\n";
        if (lowest_slope < -aor) {
          px = lowx;
          py = lowy;
          if (debug) std::cout << "    slumps to " << px << " " << py << "\n";
        } else {
          // we can stop slumping
          keep_slumping = false;
        }
      }

      // here's a final resting place
      sand(px,py) += 1;

      // update the shadow map now
      int32_t shadow_length = 0;
      int32_t shadow_height = sand(px,py);
      while (shadow_height > shadow(px,py)) {
        shadow(px,py) = shadow_height;
        shadow_height -= ashade;
        shadow_length++;
        px = (px + 1) % nx;
      }
      if (debug and shadow_length > 0) std::cout << "    casts shadow " << shadow_length << "\n";

    } // loop over test points

    //
    // write out the sand elevation matrix as a png
    //
    if (not outsand.empty()) {
      char outfile[255];
      sprintf(outfile, "%s_%04d.png", outsand.c_str(), step);

      std::cout << "\nWriting dem to " << outfile << std::endl;

      float** data = allocate_2d_array_f((int)nx, (int)ny);
      for (size_t i=0; i<nx; ++i) for (size_t j=0; j<ny; ++j) data[i][j] = sand(i,j) / (float)300.0;

      (void) write_png (outfile, (int)nx, (int)ny, FALSE, TRUE,
                        data, 0.0, 1.0, nullptr, 0.0, 1.0, nullptr, 0.0, 1.0);

      free_2d_array_f(data);
    }

  }	// loop over time steps

}
