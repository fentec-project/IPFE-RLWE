#ifndef PARAMS_H
#define PARAMS_H

static const char SIFE_Q_str[] = "2394279167817435158855681";

#define SIFE_NMODULI 3
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 44695554, 244872239156431953};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {16760833, 67043329, 2130706433};//*

// B = 4,8,16 ; l = 128,256,512
//#define SIFE_B 32
#define SIFE_B_x 4
#define SIFE_B_y 16
#define SIFE_L 785

#define SIFE_N 4096

#define SIFE_SIGMA 1

#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)

static const char SIFE_P_str[] = "50241";

//#define SIFE_SCALE_M 47655882005084197345UL //*

static const char SIFE_SCALE_M_str[]="47655882005084197345";

static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={11414086, 13786043, 483894834};	//*

#endif
