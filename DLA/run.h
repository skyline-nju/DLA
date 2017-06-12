#ifndef RUN_H
#define RUN_H
#include <vector>
#include <fstream>
#include "rand.h"
#include "comn.h"
#include "particle.h"
#include "grid.h"

void launch(Disk &p, Ran *myran);
void collision(Disk &p0, const Disk &p1, double ux, double uy, double l, double &l_hit);
bool one_step(Disk &p0, std::vector<Disk> &cluster, Ran *myran, Grid &grid);
void run(std::vector<Disk> &cluster, int nPar, Ran *myran, Grid &grid);
void cal_fractal_dim(const std::vector<Disk> &cluster);
void output(const std::vector<Disk> &cluster);
void output_xyz(const std::vector<Disk> &cluster);

#endif