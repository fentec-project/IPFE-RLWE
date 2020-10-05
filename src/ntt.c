#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "arith_rns.h"
#include "consts.h"

static const uint32_t SIFE_NTT_NINV[SIFE_NMODULI]={16756741, 67026961, 2130186241};

void poly_mul_ntt(uint32_t a[SIFE_N], uint32_t b[SIFE_N],uint32_t c[SIFE_N], uint32_t sel){

	uint32_t i;

	uint32_t a_t[SIFE_N], b_t[SIFE_N];

	for(i=0;i<SIFE_N;i++){
		a_t[i]=a[i];
		b_t[i]=b[i];
	}

	CT_forward(a_t, sel);
	CT_forward(b_t, sel);
	point_mul(a_t, b_t, c, sel);
	GS_reverse(c, sel);

}

void point_mul(uint32_t a[SIFE_N], uint32_t b[SIFE_N], uint32_t c[SIFE_N], uint32_t sel){
	
		uint64_t i;

		for(i=0;i<SIFE_N;i++){
			c[i]=mul_mod_ntt(a[i], b[i], sel);
		}
}


void CT_forward(uint32_t a[SIFE_N], uint32_t sel){


	int64_t t, m , i, j, j1, j2;
	int64_t S,U,V;

	t=SIFE_N;
	for(m=1; m<SIFE_N; m=2*m){
		t=t/2;
		for(i=0;i<m;i++){
			//printf("----\n");
			j1=2*i*t;
			j2=j1+t-1;
			//printf("Accessing twiddle factor at %ld, S : %u\n", m+i, psi[sel][m+i]);
			S=psi[sel][m+i];
			for(j=j1; j<=j2; j++){
				//printf("--\n");
				//printf("Accessing array %ld\n", j);
				//printf("Accessing array %ld\n", j+t);
				U=a[j];
				//V=reduce( (uint64_t)(a[j+t]*S), SIFE_MOD_Q_I[sel] );
				V=mul_mod_ntt( a[j+t], S, sel );
				//a[j]=reduce ( (uint64_t)(U+V), SIFE_MOD_Q_I[sel] );
				a[j]=add_mod_ntt(U, V, sel);
				//a[j+t]=reduce( (uint64_t)(U-V), SIFE_MOD_Q_I[sel] );
				a[j+t]=sub_mod_ntt(U, V, sel);
				//printf("%d, %d, %d, %d\n", a[0], a[1], a[2], a[3]);
				//print_arr(a,4);
			}
		}
	}
}


void GS_reverse(uint32_t a[SIFE_N], uint32_t sel){

	int64_t t, m, j, j1, h, i, j2;
	uint64_t S,U,V;

	//printf("In GS reverse\n");
	//printf("%d, %d, %d, %d\n", a[0], a[1], a[2], a[3]);


	t=1;
	for(m=SIFE_N; m>1; m=m/2) {
		//printf("------m=%d\n", m);
		j1=0;
		h=m/2;
		for(i=0;i<h;i++){
			//printf("----\n");
			//printf(" i: %d, h: %d\n", i, h);
			j2=j1+t-1;
			//printf("Accesing twiddle factor at %d, S=%d\n", h+i, psi_inv[sel][h+i]);
			S=psi_inv[sel][h+i];
			for(j=j1;j<=j2;j++){
				//printf("--\n");
				//printf("Accesing array %d\n", j);
				//printf("Accesing array %d\n", j+t);
				U=a[j];
				V=a[j+t];
				//a[j]=reduce( U+V, SIFE_MOD_Q_I[sel] );
				a[j]=add_mod_ntt(U, V, sel);
				//a[j+t]=reduce( (uint64_t)((U-V)*S), SIFE_MOD_Q_I[sel] );
				a[j+t]=sub_mod_ntt(U, V, sel);
				a[j+t]=mul_mod_ntt(a[j+t], S, sel);

				//printf("%d, %d, %d, %d\n", a[0], a[1], a[2], a[3]);
				//print_arr(a,4);
			}
			j1=j1+2*t;
		}
		t=2*t;
	}
		
	for(i=0;i<SIFE_N;i++){
		//a[i]=reduce( ((uint64_t)a[i]*(uint64_t)SIFE_NTT_NINV[sel]), SIFE_MOD_Q_I[sel]);	
		a[i]=mul_mod_ntt(a[i], SIFE_NTT_NINV[sel], sel);
	}
}

