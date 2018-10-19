/*
    Copyright 2002 Alfredo Braunstein, Michele Leone, Marc Mézard, 
                   Martin Weigt and Riccardo Zecchina

    This file is part of SP (Survey Propagation).

    SP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or 
    (at your option) any later version.

    SP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SP; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

//default parameters
#define   ITERATIONS 1000  // max # iterations before giving up
#define   EPSILONDEF 0.01     // convergence check
#define   NPROBE 10        // Every NPROBE fixes, run walksat (or simann)

//fixing strategy (see function fix_*)
#define   ME 0.1
#define   MEBIL 0.2
#define   MENOBIL 0.001
#define   MEZERO 0.1
#define   PARAMAGNET 0.01
#define   NDANGER 6
#define   MAXPAP 100
//#define double float

bool reached_walksat_threshold(double average_magnetization);
void run_walksat_and_exit();
void eliminate_oneclauses();
double surveyprop_iteration();
int parsecommandline(int argc, char **argv);
void randomize_eta();
void initmem();
double update_eta(int cl);
void compute_pi();
struct weightstruct compute_field(int var);
void print_fields();
void print_eta();
int run_surveyprop();
int fix(int var, int spin); 
int fix_chunk(int numInList, int quant);
void build_list(double (*)(struct weightstruct), int*, double*);
double index_biased(struct weightstruct H);
double index_best(struct weightstruct H);
double index_frozen(struct weightstruct H);
double index_para(struct weightstruct H);
double index_pap(struct weightstruct H);
double compute_sigma();
double abs_index_biased(struct weightstruct H);
double pos_index_biased(struct weightstruct H);
double neg_index_biased(struct weightstruct H);
