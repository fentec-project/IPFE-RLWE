#ifndef PARAMS_H
#define PARAMS_H

//#define SEC_LEVEL 0

#define SEC_LEVEL 1

//#define SEC_LEVEL 2


#if SEC_LEVEL==0

static const char SIFE_Q_str[] = "54453379469456060417";
#define SIFE_LOGQ_BYTE 66			//ceil(log(q)/8)
#define SIFE_NMODULI 3
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 1956705, 113253237133491840};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {12289, 8257537, 536608769};//*
#define SIFE_B_x 2
#define SIFE_B_y 2
#define SIFE_L 64
#define SIFE_N 2048
#define SIFE_SIGMA 1
#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)
static const char SIFE_P_str[] = "257";
static const char SIFE_SCALE_M_str[]="211880853966755098"; // floor(q/p)
static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={8654, 3309440, 415506400};	//*

#elif SEC_LEVEL==1

static const char SIFE_Q_str[] = "76687145727357674227351553";
#define SIFE_LOGQ_BYTE 88
#define SIFE_NMODULI 3
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 206923011, 2935204199007202395};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {16760833, 2147352577, 2130706433};//*
#define SIFE_B_x 4
#define SIFE_B_y 16
#define SIFE_L 785
#define SIFE_N 4096
#define SIFE_SIGMA 1
#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)
static const char SIFE_P_str[] = "50241";
static const char SIFE_SCALE_M_str[]="1526385735302993058007";// floor(q/p)
static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={13798054, 441557681, 1912932552};	//*

#elif SEC_LEVEL==2

static const char SIFE_Q_str[] = "1637410683940770091786553098241";
#define SIFE_LOGQ_BYTE 104
#define SIFE_NMODULI 4
static const uint64_t SIFE_CRT_CONSTS[SIFE_NMODULI]={0, 903170, 1917568766440848807, 4285461029422722560};	//*
static const uint32_t SIFE_MOD_Q_I[SIFE_NMODULI] = {114689, 1032193, 4293918721, 3221225473};//*
#define SIFE_B_x 32
#define SIFE_B_y 32
#define SIFE_L 1024
#define SIFE_N 8192
#define SIFE_SIGMA 1
#define SIFE_P (SIFE_B*SIFE_B*SIFE_L + 1)
static const char SIFE_P_str[] = "1048577";
static const char SIFE_SCALE_M_str[]="1561555025468582747653775";// floor(q/p)
static const uint64_t SIFE_SCALE_M_MOD_Q_I[SIFE_NMODULI]={17902, 926438, 3160509796, 3131397206};	//*

#endif


#endif
