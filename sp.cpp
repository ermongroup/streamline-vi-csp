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

#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "random.h"
#include "formula.h"
#include "sp.h"
#include "queue.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

//global system vars & structures
struct vstruct *v=NULL; //all per var information
struct clausestruct *clause=NULL; //all per clause information
double *magn=NULL; //spin magnetization list (for fixing ordering)
int *perm=NULL; //permutation, for update ordering
int M=0; //number of clauses
int N=0; //number of variables
int *ncl=NULL; //clause size histogram
int maxlit=0; //maximum clause size
int freespin; //number of unfixed variables
double epsilon=EPSILONDEF; //convergence criterion
int maxconn=0; //maximum connectivity

//auxiliary vars
int *list=NULL;
double *prod=NULL;

//flags & options
double percent=0;
int num_fix_per_step=1;
int generate=0;
int verbose=0;
int fields=0;
int complexity=0;
FILE *replay=NULL;
char *load=NULL;
int iterations=ITERATIONS;
int K=4;
double rho=0;
double norho=1;

int max_iter=250;
int streamlining_iter=90;
int use_streamlining=0;
int disjunction_limit=2;
int maxclauses=100000;
double ratio_bt_dec=0.0;
std::string prefix="";
int seed=0;
int proxy_freespin;
double *magn_sign=NULL;
double *pos_magn_sign=NULL;
double *neg_magn_sign=NULL;
double *free_magn_sign=NULL;

#define EPS (1.0e-16)

int main(int argc, char ** argv)
{
    parsecommandline(argc,argv);
    use_streamlining = (streamlining_iter>0);

    if (use_streamlining) {
        printf("running survey inspired streamlining\n");
    }  else {
        printf("running survey inspired decimation\n");
    }

    if (generate) {
        randomformula(K);
    } else if (load) {
        readvarformula(load);
    } else {
        fprintf(stderr, "error: you must specify some formula (-l or -n & -a)\n");
        return 0;
    }

    //define num_fix_per_step if based on percent
    if(percent) {
        num_fix_per_step=N*percent/100.0;
    }
    maxclauses = M + 100*N;
    initmem();
    // Simplify formula 
    eliminate_oneclauses();

    std::string NEW_FORMULA = prefix + FORMULA;
    writeformula(fopen(NEW_FORMULA.c_str(),"w+"));
    randomize_eta();
    std::string STREAMLINE_FORMULA = prefix + STREAMLINEFORMULA;
    
    for (int iterNum=0; iterNum<max_iter; iterNum++) {    
        if (!run_surveyprop()) {
            exit(-1);
        }
        int oldfreespin = freespin;

        int numInList;
        double average_magnetization;
        build_list(&index_biased, &numInList, &average_magnetization);
        if (reached_walksat_threshold(average_magnetization) or iterNum == max_iter-1) {
            run_walksat_and_exit();
        }

        if (use_streamlining && 
                iterNum > streamlining_iter) {
            use_streamlining = 0;
            printf("streamlining finished...decimating now\n");
            writeformula(fopen(STREAMLINE_FORMULA.c_str(),"w+"));
        }
        
        fix_chunk(numInList, num_fix_per_step);

        if(verbose){
            printf("fixed %i biased var%s (+%i ucp)\n",num_fix_per_step,num_fix_per_step>1?"s":"",oldfreespin-freespin-num_fix_per_step);
            print_stats();
        }
    }
    return -1;
}

bool reached_walksat_threshold(double average_magnetization) {
    return average_magnetization < PARAMAGNET;
}

void run_walksat_and_exit() {
    printf("\nparamagnetic state, small s case\n");
    printf("sub-formula has:\n");
    print_stats();
    std::string NEW_SPSOL = prefix + SPSOL;
    writesolution(fopen(NEW_SPSOL.c_str(),"w+"));
    solvesubsystem();
    exit(0);
}

void eliminate_oneclauses() {
    int c,var;
    for(c=0; c<M; c++) { 
        if(clause[c].lits==1) {
            var=clause[c].literal[0].var;
            if(v[var].spin==0) {
                printf("eliminating var %i (in 1-clause)\n",var);
                fix(var,clause[c].literal[0].bar?-1:1);
            }
        }
    }
}

/************************************************************
 ************************************************************
 * VARIABLE SELECTION
 ************************************************************
 ************************************************************/

inline struct weightstruct normalize(struct weightstruct H)
    //normalize a triplet
{
    double norm;
    norm=H.m+H.z+H.p;
    H.m/=norm;
    H.p/=norm;
    H.z/=norm;
    return H;
}

int order(void const *a, void const *b)
    //order relation for qsort, uses ranking in magn[]
{
    double aux;
    aux=magn[*((int *)b)]-magn[*((int *)a)];
    return aux<0?-1:(aux>0?1:0);
}

double index_biased(struct weightstruct H)
    //most biased ranking
{
    return fabs(H.p-H.m);
}
double index_pap(struct weightstruct H)
    //most biased with some fuss
{
    return fabs(H.p-H.m)+randreal()*0.1;
}

double index_para(struct weightstruct H)
    //least paramagnetic ranking
{
    return H.z;
}

double index_frozen(struct weightstruct H)
    //most paramagnetic ranking
{
    return -H.z;
}

double index_best(struct weightstruct H)
    //min(H.p,H.m) ranking
{
    return 1.-(H.p>H.m?H.m:H.p);
}


// added
// Following our discussion, taking H.m as positive and H.p as negative
double abs_index_biased(struct weightstruct H)
    //abs biased ranking
{
    return H.m-H.p;
}


double pos_index_biased(struct weightstruct H)
    //abs biased ranking
{
    return H.m+H.z/2;
}


double neg_index_biased(struct weightstruct H)
    //abs biased ranking
{
    return H.p+H.z/2;
}


void build_list(
        double (* index)(struct weightstruct),
        int *numInList,
        double *averageMagnetization
        )
    //build an ordered list with order *index which is one of index_?
{
    int var;
    struct weightstruct H;
    double summag;
    double maxmag;
    summag=0;
    *numInList=0;
    for(var=1; var<=N; var++) {
        if(v[var].spin==0 && v[var].proxy_spin==0) {
            H=compute_field(var);
            list[(*numInList)++]=var;
            magn[var]=(*index)(H);
            maxmag=H.p>H.m?H.p:H.m;

            summag+=maxmag;
            }
            if(v[var].spin==0){
            // added
            magn_sign[var]=abs_index_biased(H);
            pos_magn_sign[var]=pos_index_biased(H);
            neg_magn_sign[var]=neg_index_biased(H);
            free_magn_sign[var]=index_para(H);
            }
    }
    qsort(list, (*numInList), sizeof(int), &order);
    *averageMagnetization = summag / (*numInList);
}

void filterCandidateList(int numInList, std::vector<int>& out) {
    out.clear();
    for (int i = 0; i < numInList; i++) {
        int var = list[i];
        if (var <=0 || var > N) continue;
        if (v[var].spin != 0) continue;
        if(use_streamlining && (v[var].artificialclauses>disjunction_limit)) continue;
        out.push_back(var);
    }
}

int fix_chunk(int numInList, int numToFix) {
    std::vector<int> candidates;
    filterCandidateList(numInList, candidates);

    std::cout << numToFix << " " << candidates.size() << std::endl;

    int last_var = std::min(
            (size_t) (use_streamlining ? numToFix * 2 : numToFix), 
            candidates.size());

    for (int i = 0; i < last_var; i++) {
        int var1 = candidates[i];
        struct weightstruct H1 = compute_field(var1);
        if(verbose>1) {
            printf("H1(%i)={%f,%f,%f} ---> %s\n",var1,
                    H1.p,H1.z,H1.m,H1.p>H1.m?"-":"+");
        }

        if (use_streamlining) {
            int var2_ind = last_var - 1 - i;
            if (var2_ind <= i) break;
            int var2 = candidates[var2_ind];

            struct weightstruct H2 = compute_field(var2);
            double m1 = fabs(H1.p - H1.m);
            double m2 = fabs(H2.p - H2.m);

            if(verbose>1) {
                printf("H2(%i)={%f,%f,%f} ---> %s\n",var2,
                        H2.p,H2.z,H2.m,H2.p>H2.m?"-":"+");
            }

            if (verbose > 1) {
                    std::cout << "OR constraining " 
                        << var1 << " (m=" << m1 << ") and " 
                        << var1 << " (m=" << m2 << ")" << std::endl;
            }
            fixTwo(var1, var2, H1.p>H1.m?-1:1, H2.p>H2.m?-1:1);
        } else {
            fix(var1,H1.p>H1.m?-1:1);
        }

        // Once we've fixed the last needed clause, break.
        if (numToFix-- == 1) break;
    } 
    return numToFix;
}

/************************************************************
 ************************************************************
 * SURVEY PROPAGATION
 ************************************************************
 ************************************************************/

void randomize_eta()
    //pick initial random values
{
    int i,j;
    for(i=0; i<M; i++) {
        for(j=0; j<clause[i].lits;j++) {
            clause[i].literal[j].eta=randreal();
        }
    }
}

void initmem()
    //allocate mem (can be called more than once)
{
    delete[] perm;
    delete[] list;
    delete[] magn;
    delete[] prod;
    delete[] magn_sign; //added
    delete[] pos_magn_sign; //added
    delete[] neg_magn_sign; //added
    delete[] free_magn_sign; //added

    perm=new int[maxclauses];
    list=new int[N+1];
    magn=new double[N+1];
    prod=new double[maxlit];
    magn_sign = new double[N+1]; //added
    pos_magn_sign = new double[N+1]; //added
    neg_magn_sign = new double[N+1]; //added
    free_magn_sign = new double[N+1]; //added

    queue_init(M);
    if(!perm||!list||!magn ||!prod||!magn_sign || !pos_magn_sign || !neg_magn_sign || !free_magn_sign){
        fprintf(stderr, "Not enough memory for internal structures\n");
        exit(-1);
    }
}

void compute_pi()
    //compute pi products of all vars from scratch
{
    int i,var;
    struct clauselist *cl;
    struct literalstruct *l;
    for(var=1; var<=N; ++var) if(v[var].spin==0) {
        v[var].pi.p=1;
        v[var].pi.m=1;
        v[var].pi.pzero=0;
        v[var].pi.mzero=0;
        for(i=0,cl=v[var].clauselist;i<v[var].clauses; i++,cl++) if(cl->clause->type) {
            l=cl->clause->literal+cl->lit;
            if(l->bar) {
                if(1-l->eta>EPS) {
                    v[var].pi.p*=1-l->eta;
                } else {
                    v[var].pi.pzero++;
                }
            } else {
                if(1-l->eta>EPS) {
                    v[var].pi.m*=1-l->eta;
                } else {
                    v[var].pi.mzero++;
                }
            }
        }
    }
}

double update_eta(int cl)
    // updates all eta's and pi's in clause cl 
{
    int i;
    struct clausestruct *c;
    struct literalstruct *l;
    struct pistruct *pi;
    double eps;
    double neweta;
    double allprod=1;
    double wt,wn;
    int zeroes=0;
    double m,p;

    c=&(clause[cl]);
    for(i=0,l=c->literal; i<c->lits; i++,l++) if(v[l->var].spin==0) {
        pi=&(v[l->var].pi);
        if(l->bar) {
            m=pi->mzero?0:pi->m;
            if(pi->pzero==0) {
                p=pi->p/(1-l->eta);
            } else if (pi->pzero==1 && 1-l->eta<EPS) {
                p=pi->p;
            } else {
                p=0;
            }
            wn=p*(1-m*norho);
            wt=m;
        } else { 
            p=pi->pzero?0:pi->p;
            if(pi->mzero==0) {
                m=pi->m/(1-l->eta);
            } else if (pi->mzero==1 && 1-l->eta<EPS) {
                m=pi->m;
            } else {
                m=0;
            }
            wn=m*(1-p*norho);
            wt=p;
        }
        prod[i]=wn/(wn+wt);
        if(prod[i]<EPS) {
            if(++zeroes == 2)
                break;
        } else {
            allprod*=prod[i];
        }
    }
    eps=0;
    for(i=0,l=c->literal; i<c->lits; i++,l++) if(v[l->var].spin==0) {
        if(!zeroes){
            neweta=allprod/prod[i];
        } else if (zeroes==1 && prod[i]<EPS) {
            neweta=allprod;
        } else {
            neweta=0;
        }

        pi=&(v[l->var].pi);
        if(l->bar) {
            if(1-l->eta>EPS) {
                if(1-neweta>EPS) {
                    pi->p*=(1-neweta)/(1-l->eta);
                } else {
                    pi->p/=(1-l->eta);
                    pi->pzero++;
                } 
            } else {
                if(1-neweta>EPS) {
                    pi->p*=(1-neweta);
                    pi->pzero--;
                } 
            }
        } else {
            if(1-l->eta>EPS) {
                if(1-neweta>EPS) {
                    pi->m*=(1-neweta)/(1-l->eta);
                } else {
                    pi->m/=(1-l->eta);
                    pi->mzero++;
                } 
            } else {
                if(1-neweta>EPS) {
                    pi->m*=(1-neweta);
                    pi->mzero--;
                } 
            }
        }
        if(eps<fabs(l->eta-neweta))
            eps=fabs(fabs(l->eta-neweta));
        l->eta=neweta;
    }
    return eps;
}

struct weightstruct compute_field(int var)
    //compute H field of the variable var
{
    struct weightstruct H;
    double p,m;
    p=v[var].pi.pzero?0:v[var].pi.p;
    m=v[var].pi.mzero?0:v[var].pi.m;
    H.z=p*m;
    H.p=m-H.z;
    H.m=p-H.z;
    return normalize(H);
}

void print_fields()
    //print all H (non-cavity) fields
{
    int var;
    struct weightstruct H;
    compute_pi();
    for(var=1; var<=N; var++) if(v[var].spin==0) {
        H=compute_field(var);
        printf("#H(%i)={%f,%f,%f}\n",var,H.p,H.z,H.m);
    }
}

void print_eta()
    //print all etas
{
    int c,l;
    for(c=0; c<M; c++)  if(clause[c].type){
        for(l=0;l<clause[c].lits;l++) if(!v[clause[c].literal[l].var].spin){
            printf("#eta(%i,%i)=%f\n",c,l,clause[c].literal[l].eta);
        }
    }
}


double compute_sigma()
    //compute complexity
{
    int var,cl,i,n;
    double allprod,allprod1,wt,wn;
    double sigmabonds=0, sigmasites=0, sigma=0;
    struct literalstruct *l;
    struct clausestruct *c;
    struct pistruct *pi;
    struct clauselist *cli;
    sigmabonds=0;
    for(cl=0, c=clause; cl<M; c++, cl++) if(c->type) {
        allprod=1;
        allprod1=1;
        for(i=0,l=c->literal; i<c->lits; i++,l++) if(v[l->var].spin==0) {
            pi=&(v[l->var].pi);
#ifdef FAST_ITERATION
            if(l->bar) {
                wn=pi->p/(1-l->eta);
                wt=pi->m;
            } else {
                wn=pi->m/(1-l->eta);
                wt=pi->p;
            }
#else
            if(l->bar) {
                if(1-l->eta>EPS) {
                    wn=pi->p/(1-l->eta);
                } else if (pi->pzero==1) {
                    wn=pi->p;
                } else {
                    wn=0;
                }
                wt=pi->mzero?0:pi->m;
            } else { 
                if(1-l->eta>EPS) {
                    wn=pi->m/(1-l->eta);
                } else if (pi->mzero==1) {
                    wn=pi->m;
                } else {
                    wn=0;
                }
                wt=pi->pzero?0:pi->p;
            }
#endif
            allprod*=wn*(1-wt);
            allprod1*=wt+wn-wt*wn;
        }
        sigmabonds+=log(allprod1-allprod);
    }
    sigmasites=0;
    for(var=1;var<=N;var++) if(v[var].spin==0){
        n=0;
        for(cli=v[var].clauselist,cl=0; cl<v[var].clauses; cli++,cl++) {
            if(cli->clause->type)
                n++;
        }
#ifdef FAST_ITERATION
        wt=v[var].pi.p;
        wn=v[var].pi.m;
#else
        wt=v[var].pi.pzero?0:v[var].pi.p;
        wn=v[var].pi.mzero?0:v[var].pi.m;
#endif
        sigmasites+=log(wt+wn-wt*wn)*(n-1);
    }
    sigma=(sigmabonds-sigmasites)/freespin;
    printf("bonds=%f sites=%f sigma=%f\n", sigmabonds,sigmasites,sigma);
    return sigma;
}

int run_surveyprop() {
    /* Run surveyprop iterations until convergence. */
    double eps=0;
    int iter=0,cl,quant;
    compute_pi();

    for(cl=0,quant=0; cl<M; cl++) if(clause[cl].type) {
        perm[quant++]=cl;
    }

    do {
        eps=surveyprop_iteration();
        fflush(stderr);
    } while (eps>epsilon && iter++<iterations);
    if(eps<=epsilon) {
        return 1;
    } else {
        printf("SP did NOT converge for seed:%d t:%d\n", seed, streamlining_iter);
        writeformula(fopen(NOCONVERGENCE,"w+"));
        return 0;
    }

}

double surveyprop_iteration() {
    int cl,vi=0,quant,i;

    double eps,maxeps;
    eps=0;
    maxeps=0;
    for(quant=M-ncl[0];quant;--quant) {
        cl=perm[i=randint(quant)];
        perm[i]=perm[quant-1];
        perm[quant-1]=cl;
        eps=update_eta(cl);
        if(eps>epsilon) {
            vi++;
        }
        if(eps>maxeps) {
            maxeps=eps;
        }
    }
    return maxeps;
}

int parsecommandline(int argc, char **argv)
    //get command line parameters
{   
    double alpha=0;
    char c;
    generate=0;
    usrand(1);
    extern char* optarg;
    while((c=getopt(argc, argv,
                    "R:k:cN:M:r:n:m:a:s:hf:v%:e:l:dFi:Q:t:p:"))!=-1) {
        switch (c) {
            case 'R':
                rho=atof(optarg);
                norho=1-rho;
                break;
            case 'l':
                load=optarg;
                break;
            case 'r':
                replay=fopen(optarg,"w+");
                break;
            case 'N':
            case 'n':
                N=atoi(optarg);
                break;
            case 'M':
            case 'm':
                M=atoi(optarg);
                break;
            case 'a':
                alpha=atof(optarg);
                break;
            case 'e':
                epsilon=atof(optarg);
                break;
            case 't':
                streamlining_iter=atoi(optarg);
                break;
            case 's':
                seed=atoi(optarg);
                usrand(atoi(optarg));
                srand(atoi(optarg));
                break;
            case 'k':
            case 'K':
                K=atoi(optarg);
                break;
            case 'v':
                verbose++;
                break;
            case 'F':
                fields=1;
                break;
            case 'f':
                num_fix_per_step=atoi(optarg);
                break;
            case '%':
                percent=atof(optarg);
                break;
            case 'i':
                iterations=atoi(optarg);
                break;
            case 'Q':  
                ratio_bt_dec=(double) atof(optarg);
                break;
            case 'c':
                complexity=1;
                break;
            case 'd':
                disjunction_limit=atoi(optarg);
                break;
            case 'p':
                prefix=optarg;
                break;
            case 'h':
            default:
                fprintf(stderr, "%s [options]\n"
                        "\n"
                        "  formula\n"
                        "\t-n <numvars>\n"
                        "\t-m <numclauses>\n"
                        "\t-a <alpha>\n"
                        "\t-R <rho>\t modified dynamics (0=sp, 1=bp)\n"
                        "\t\t\t (real values inbetween may make sense)\n"
                        "\t-l <filename>\t reads formula from file\n"
                        "  solving\n"
                        "\t-f <fixes>\t per step\n"
                        "\t-%% <fixes>\t per step (%%)\n"
                        "\t-e <error>\t criterion for convergence\n"
                        "\t-z \t\t use lazy convergence instead of sequential\n"
                        "\t-i <iter>\t maximum number of iterations for survey prop until convergence\n"
                        "  stats\n"
                        "\t-p <filename>\t output a magneticity plot\n"
                        "\t-r <filename>\t replay file\n"
                        "\t-c \t\t computes complexity\n"
                        "\t-F \t\t print fields\n"
                        "\t-v \t\t increase verbosity\n"
                        "  misc\n"
                        "\t-s <seed>\t (0=use time, default=1)\n"
                        "\t-t <iter>\t streamlining iterations (0: decimation only, >0: streamlining for t iterations before decimation)\n"
                        "\t-d <disjunctionlimit>\t limit on the number of extra disjunctions per variable\n"
                        "\t-h\t\t this help\n",argv[0]);
                exit(-1);
        }
    }
    if(load && !N && !M && !alpha) {
        generate=0;
    } else if(N && alpha && !M) {
        M=N*alpha;
        generate=1;
    } else if(M && alpha && !N) {
        N=M/alpha;
        generate=1;
    } else if(M && N && alpha==0) {
        generate=1;
    } else {
        fprintf(stderr, "error: you have to specify exactly TWO of -n,-m and -a, or -l FILE (and then a formula is read from FILE)\n");
        exit(-1);
    }
    return 0;
}


