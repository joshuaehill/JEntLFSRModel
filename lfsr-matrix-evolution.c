#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>

/* This program enumerates a range of expected raw values, and runs the LFSR forward for each,
 * saving the outputs to support statistical testing of them.
 */

#define OUTPUTBLOCKS UINT64_C(125000) // Each block is 8 bytes, so this represents 1 million bytes

//Narrow results (ARMish)
/*
#define TESTBITWIDTH UINT64_C(3)
#define TESTOFFSET UINT64_C(120)
*/

//Wide results (x86-with-TSCish)
#define TESTBITWIDTH UINT64_C(8)
#define TESTOFFSET UINT64_C(8645)

//Some conventions:
//Column vectors are stored as uint64_t values
//The standard basis is e_1 = 0x0000000000000001, e_2 = 0x0000000000000002, e_3 = 0x0000000000000004, etc.

// This is the full LFSR process where the input value is shifted into the internal state.
// It is derived from JEnt v2.2.0 jent_lfsr_time() function.
// The changes here are to isolate the LFSR processing from the other JEnt functionality
// and to harmonize the variable naming with the LFSR conditioning analysis.
static uint64_t lfsr_process(uint64_t s, uint64_t x) {

	//Shifting in the input (x) low order to high order
	for (uint64_t i = 1; 64 >= i; i++) {
		//This operation does a (64-i) left bit shift;
		//this puts a shifted copy of x into tmp which
		//makes bit place (i-1) of x the msb of tmp.
		//(The loop left-shifts 63 then 62 then 61 ... and finally 0 bits)
		uint64_t tmp = x << (64 - i);
		//shift the 64-bit value tmp right by 63 bits, 
		//retaining only the msb of tmp 
		//(which is the (i-1)th bit of x) 
		//and positioning it in the lsb of tmp.
		tmp = tmp >> 63;

		//Now tmp contains only the (i-1)th bit of x in its lsb.

		/*
		* Fibonacci LFSR with feedback polynomial of
		*  z^64 + z^61 + z^56 + z^31 + z^28 + z^23 + 1 
		* which is primitive according to
		*   http://poincare.matf.bg.ac.rs/~ezivkovm/publications/primpol1.pdf
		* (the shift values are the polynomial powers minus one
		* due to counting bits from 0 to 63). As the current
		* position is always the lsb, the polynomial only needs
		* to shift data in from the left without wrap.
		*/
		tmp ^= ((s >> 63UL) & 1UL);
		tmp ^= ((s >> 60UL) & 1UL);
		tmp ^= ((s >> 55UL) & 1UL);
		tmp ^= ((s >> 30UL) & 1UL);
		tmp ^= ((s >> 27UL) & 1UL);
		tmp ^= ((s >> 22UL) & 1UL);
		s <<= 1UL;
		s ^= tmp;
	}

	return s;
}

//Some basic matrix math functions.
//We're trying for clarity over speed here.
// Return the result of the matrix multiplication M v as a column matrix encoded as a uint64_t
// We can calculate each term of M v = y as y _i = Sum_{k=1}^64 M_{i,k} v_k
// Here, we are operating in a binary ring, so addition is XOR
// Remember that C indexes start at 0, whereas the standard way of writing matrix column/row
// indices is to start at 1.
static uint64_t applyMatrix(uint8_t M[][64], uint64_t v) {
	uint64_t output = 0;

	for(uint8_t i=0; i<64; i++) {
		uint64_t curres = 0;
		//Calculate the result for row i+1 of the result
		//Recall that C is 0 indexed, whereas matrix indexes are 1-indexed.
		//(e.g., the i=0 case is for row 1 of the matrix)
		for(uint8_t k=0; k<64; k++) {
			uint8_t vk = ((v&(1UL<<k))==0UL)?0:1;
			curres = curres ^ (M[i][k] * vk);
		}
		assert((curres & 0xFFFFFFFFFFFFFFFE)==0UL);
		//We now have calculated y_i; it resides in the lsb of curres.
		output = output ^ (curres << i);
	}
	return output;
}

int main(void) {
	uint8_t A[64][64];
	uint8_t B[64][64];

	//Calculate the linear maps by applying the functions to the standard basis.
	//A
	for(uint64_t j=0; j<64; j++) {
		uint64_t result = lfsr_process(1UL<<j, 0UL);
		for(uint64_t i=0; i<64; i++) A[i][j] = ((result&(1UL<<i))==0UL)?0:1;
	}

	//B
	for(uint64_t j=0; j<64; j++) {
		uint64_t result = lfsr_process(0UL, 1UL<<j);
		for(uint64_t i=0; i<64; i++) B[i][j] = ((result&(1UL<<i))==0UL)?0:1;
	}

	for (uint64_t x=TESTOFFSET; x < TESTOFFSET+(UINT64_C(1)<<TESTBITWIDTH); x++) {
		char curfilename[256];
		FILE *fp;
		uint64_t curvalue=applyMatrix(B, x);

		printf("Testing initial value %" PRIu64 "\n", x);

		snprintf(curfilename, sizeof(curfilename), "out-%" PRIu64 ".bin", x);
		
		if((fp=fopen(curfilename, "wb"))==NULL) {
			perror("Can't open file.");
			exit(-1);
		}

		//Write out the initial value B x (a.k.a. A^0 B x)
		if(fwrite(&curvalue, sizeof(curvalue), 1, fp)!=1) {
			perror("Can't write to file.");
			exit(-1);
		}

		for(uint64_t j = 1; j < OUTPUTBLOCKS; j++) {
			curvalue=applyMatrix(A,curvalue);
			//write out A^j B x
			if(fwrite(&curvalue, sizeof(curvalue), 1, fp)!=1) {
				perror("Can't write to file.");
				exit(-1);
			}
		}
		fclose(fp);
	}

	return 0;
}
