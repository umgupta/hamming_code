#include<stdio.h>
#include<sys/time.h>
#include<arm_neon.h>
#include<stdint.h>
#include<stdlib.h>
#include<string.h>
#define SIZE 60000
#define OFFSET 0
#define USAGE_ARG1 "./a.out arg1 arg2 \n arg1 -\tgive 1 to run a 16-bit implementation, \n\tgive 2 to run a 64-bit implementation \n "
#define USAGE_ARG2 "arg2-\tgive 1 to align memory to 64 byte \n\tgive 2 to run with unaligned memory (unaligned to 64 bytes)\n"
//you have to pass offset and size in units of byte/ size is a multiple of 8 and offset is a multiple of 2.
static inline int compute_ham_similarity_64(unsigned short* , unsigned short*,
                                int);
static inline int compute_ham_similarity_16(unsigned short* , unsigned short*,
                                int, int);

static inline int setbits(uint16x8_t );

int main(int argc, char *argv[]){
	
if (argc!=3){
		printf (USAGE_ARG1);
		printf (USAGE_ARG2);
		exit (0);
	}
	
	int i,count,size=SIZE,offset=OFFSET;
	printf("size in bytes is %d and offset is %d\n",size,offset);
	unsigned short* reference;
	unsigned short* circ;
    struct timeval start, end;

    if(0==strcmp(argv[2],"2")){
		printf("running unaligned \n");		
		reference=malloc(size);
		circ=malloc(size);        
	}else if (0==strcmp(argv[2],"1")){
		printf("running aligned to 64 bytes\n");
		reference=malloc(size+64);
		reference=(unsigned short*)(((int)reference+64)&(~0x3F));
			
		circ=malloc(size+64);
		circ=(unsigned short*)(((int)circ+64)&(~0x3F));
	}else{
		printf("arg2 should be 1 or 2 only\n");
		exit(0);
	}


	for (i=0;i<size/2;i++){
                reference[i]=1;//*random();
                circ[i]=2;//*random();
        }
	if(0==strcmp(argv[1],"1")){
		printf("running 16 bit version\n");
		gettimeofday(&start,NULL);
       	count=compute_ham_similarity_16(reference,circ,offset/2,size/2);
   	 	gettimeofday(&end,NULL);
	}else if(0==strcmp(argv[1],"2")){
		printf("running 64 bit version\n");
		gettimeofday(&start,NULL);
        count=compute_ham_similarity_64(reference,circ,size);
        gettimeofday(&end,NULL);
	}else{
		printf("-----------------\n");
		exit(0);
	}

        printf("%ld micro-seconds elapsed\n", end.tv_sec*1000000+(end.tv_usec)-start.tv_sec*1000000-(start.tv_usec));
	printf("Hamming value: %d\n",count);
}

static inline int setbits(uint16x8_t d){
        
	return (int) (vgetq_lane_u16 (d,0)+ vgetq_lane_u16 (d,1)+ vgetq_lane_u16 (d,2)+
		      vgetq_lane_u16 (d,3)+ vgetq_lane_u16 (d,4)+ vgetq_lane_u16 (d,5)+
		      vgetq_lane_u16 (d,6)+ vgetq_lane_u16 (d,7) );
}

//Note: it takes size and offset in units of byte
static inline int compute_ham_similarity_64(unsigned short* ref, unsigned short* circ_array, 
                                int size){
	
    const uint8_t*         	ref_c=(uint8_t*) ref;
    const uint8_t*         	circ_c=(uint8_t*) circ_array;
    register uint8x16_t     a,b;
    register uint8x16_t     c,d,temp;
    register uint16x8_t     acc;
    register uint           i=0,count=0;
	int 					j=0;
	int 					shift=size&0xF;
	for(i=0;i<=size-16; i+=16){
		j++;
		a=vld1q_u8(&ref_c[i]);
		b=vld1q_u8(&circ_c[i]);
		c=veorq_u8(a,b);
		acc=vaddq_u16(acc,vpaddlq_u8(vcntq_u8(c)));
	}	
	count=setbits(acc);
	a=vld1q_u8(&ref_c[i]);
	b=vld1q_u8(&circ_c[i]);
	c=veorq_u8(a,b);
	c=vcntq_u8(c);
	
	for(i=0;i<shift;i++){
		count=count+vgetq_lane_u8 (c,i);
	}
    return size*8-count;
}

static const unsigned char BitsSetTable256[256] =
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
   B6(0), B6(1), B6(1), B6(2)
};
//Note: This takes size and offset in units of half word
static inline int compute_ham_similarity_16(unsigned short *ref, unsigned short *circ_array,
                       int off, int  size){
        int              total_bits_diff=0;
        int              i = 0;
        unsigned short   val = 0;
        int              c = 0;
        unsigned char    *p = (unsigned char *) &val;
        for (i=0; i<size-off; i++){
               	val = ref[i] ^ circ_array[i+off];
               	c = BitsSetTable256[p[0]] + BitsSetTable256[p[1]] ;
               	total_bits_diff += c;
	}
        for (i=size-off; i<size; i++) {
               	val = ref[i] ^ circ_array[i-(size-off)];
               	c = BitsSetTable256[p[0]] + BitsSetTable256[p[1]] ;
               	total_bits_diff += c;
        }
        return ((size)*2*8-total_bits_diff);
}

