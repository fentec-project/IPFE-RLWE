#ifndef PARAMS_H
#define PARAMS_H

// log q = 48,54,61 for current parameter sets
#define MIFE_Q1 15485863
#define MIFE_Q2 16777213
// Q = Q1*Q2 = 259809622039819
//#define MIFE_Q (MIFE_Q1*MIFE_Q2)
#define MIFE_Q 259809622039819

//C2 is a constante used for Garner's algorithm for CRT^(-1)
#define MIFE_C2 9157109

//inverses of y_i = Q/q_i
#define MIFE_Z1 7033581
#define MIFE_Z2 9157109
//coefficients for a_i = MIFE_Q2 * MIFE_Z1 mod MIFE_Q
#define MIFE_CRT_X1 118003886589753
#define MIFE_CRT_X2 141805735450067

#define MIFE_NMODULI 2

static const uint32_t MIFE_MOD_Q_I[MIFE_NMODULI] = {MIFE_Q1, MIFE_Q2};
static const uint32_t MIFE_CRT_Y_I[MIFE_NMODULI] = {MIFE_Q2, MIFE_Q1};
//static const uint32_t MIFE_CRT_Z_I[MIFE_NMODULI] = {MIFE_Z1, MIFE_Z2};
static const uint64_t MIFE_CRT_REVERSE_I[MIFE_NMODULI] = {MIFE_CRT_X1, MIFE_CRT_X2};

// B = 4,8,16 ; l = 128,256,512
#define MIFE_BBITS 2
#define MIFE_B (1 << MIFE_BBITS)
#define MIFE_L 128

#define MIFE_K 11
#define MIFE_N (1 << MIFE_K)

//#define MIFE_SIGMA 2*MIFE_K
#define MIFE_SIGMA 1
//epsilon = 2^(-80)
#define MIFE_EPSILON_EXP -80
//p > 2*B²*l
#define MIFE_P (2*MIFE_B*MIFE_B*MIFE_L + 1)
//floor(q/p)
#define MIFE_SCALE_M (MIFE_Q/MIFE_P)

static const uint64_t MIFE_SCALE_M_MOD_Q_I[MIFE_NMODULI] = {(MIFE_SCALE_M % MIFE_Q1), (MIFE_SCALE_M % MIFE_Q2)};

//sigma1 = sqrt(2l)*B*sigma
//sigma2 = sqrt(2)*sigma
//sigma3 = sqrt(2(n*l*sigma2²*log(1/epsilon)+1)*sigma1
#define MIFE_SIGMA1 10
#define MIFE_SIGMA2 10
#define MIFE_SIGMA3 10

#endif