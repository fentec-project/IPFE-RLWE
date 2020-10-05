#ifndef PARAMS_H
#define PARAMS_H

//#define SEC_LEVEL 0

//#define SEC_LEVEL 1

#define SEC_LEVEL 2


#if SEC_LEVEL==0

static const char SIFE_Q_str[] = "15131549376305999873";
#define SIFE_NMODULI 3
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 70580, 1461855532664819506};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {61441, 114689, 2147352577};//*
#define SIFE_B_x 2
#define SIFE_B_y 2
#define SIFE_L 64
#define SIFE_N 2048
#define SIFE_SIGMA 1
#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)
static const char SIFE_P_str[] = "257";
static const char SIFE_SCALE_M_str[]="58877624032319065";
static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={22711, 110672, 1487271434};	//*

#elif SEC_LEVEL==1

static const char SIFE_Q_str[] = "2394279167817435158855681";
#define SIFE_NMODULI 3
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 44695554, 244872239156431953};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {16760833, 67043329, 2130706433};//*
#define SIFE_B_x 4
#define SIFE_B_y 16
#define SIFE_L 785
#define SIFE_N 4096
#define SIFE_SIGMA 1
#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)
static const char SIFE_P_str[] = "50241";
static const char SIFE_SCALE_M_str[]="47655882005084197345";
static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={11414086, 13786043, 483894834};	//*

#elif SEC_LEVEL==2

static const char SIFE_Q_str[] = "19796162090054631350254043137";
#define SIFE_NMODULI 3
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 2146953901, 8460412124};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {2147352577, 2146959361 , 4293918721};//*
#define SIFE_B_x 32
#define SIFE_B_y 32
#define SIFE_L 1024
#define SIFE_N 8192
#define SIFE_SIGMA 1
#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)
static const char SIFE_P_str[] = "1048577";
static const char SIFE_SCALE_M_str[]="18879073344212805879066";
static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={2126372117, 1592386324, 1319825431};	//*

#endif


#endif
