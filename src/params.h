#ifndef PARAMS_H
#define PARAMS_H

static const char MIFE_Q_str[] = "2394279167817435158855681";

#define MIFE_NMODULI 3
static const uint64_t MIFE_CRT_CONSTS[MIFE_NMODULI]={0, 44695554, 244872239156431953};	//*
static const uint32_t MIFE_MOD_Q_I[MIFE_NMODULI] = {16760833, 67043329, 2130706433};//*

// B = 4,8,16 ; l = 128,256,512
//#define MIFE_B 32
#define MIFE_B_x 4
#define MIFE_B_y 16
#define MIFE_L 785

#define MIFE_N 4096

#define MIFE_SIGMA 1

#define MIFE_P (MIFE_B*MIFE_B*MIFE_L + 1)

static const char MIFE_P_str[] = "50241";

//#define MIFE_SCALE_M 47655882005084197345UL //*

static const char MIFE_SCALE_M_str[]="47655882005084197345";

static const uint64_t MIFE_SCALE_M_MOD_Q_I[MIFE_NMODULI]={11414086, 13786043, 483894834};	//*


#define MIFE_SIGMA1 10
#define MIFE_SIGMA2 10
#define MIFE_SIGMA3 10

#endif
