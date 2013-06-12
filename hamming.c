#include<stdio.h>
#include<sys/time.h>
#include<arm_neon.h>
#include<stdint.h>
#include<stdlib.h>
#include<string.h>
#define USAGE "Give 1 for 16 bit implementation, 2 for 64-bit implementation"
//you have to pass offset and size in units of byte/ size is a multiple of 8 and offset is a multiple of 2.
static inline int compute_ham_similarity_64(unsigned short* , unsigned short*,
                                int, int);
static inline int compute_ham_similarity_16(unsigned short* , unsigned short*,
                                int, int);

static inline int setbits(uint64x1_t );

int main(int argc, char *argv[]){
	
	if (argc!=2){
	printf (USAGE);
	exit (0);
	}
	int i,count,size=600,offset=100;
	printf("size in bytes is %d and offset is %d\n",size,offset);

        struct timeval start, end;
        unsigned short reference[size/2], circ[size/2];

	for (i=0;i<size/2;i++){
                reference[i]=random();
                circ[(i+offset/2)%(size/2)]=random();
        }
	if(0==strcmp(argv[1],"1")){
		printf("running 16 bit version\n");
		gettimeofday(&start,NULL);
        	count=compute_ham_similarity_16(reference,circ,offset/2,size/2);
       	 	gettimeofday(&end,NULL);
	}else if(0==strcmp(argv[1],"2")){
		printf("running 64 bit version\n");
		gettimeofday(&start,NULL);
                count=compute_ham_similarity_64(reference,circ,offset,size);
                gettimeofday(&end,NULL);
	}else{
		printf("-----------------\n");
		exit(0);
	}

        printf("%ld micro-seconds elapsed\n", end.tv_sec*1000000+(end.tv_usec)-start.tv_sec*1000000-(start.tv_usec));
	printf("Hamming value: %d\n",count);
}

static inline int setbits(uint64x1_t a){
        
	register uint8x8_t d= vcnt_u8(vreinterpret_u8_u64 (a));
	return (int) (vget_lane_u8 (d,0)+ vget_lane_u8 (d,1)+ vget_lane_u8 (d,2)+
		      vget_lane_u8 (d,3)+ vget_lane_u8 (d,4)+ vget_lane_u8 (d,5)+
		      vget_lane_u8 (d,6)+ vget_lane_u8 (d,7) );
}

//Note: it takes size and offset in units of byte
static inline int compute_ham_similarity_64(unsigned short* ref, unsigned short* circ_array, 
                                int off, int size){

	const uint64_t* 	ref_w=(uint64_t*) ref;
        const uint64_t* 	circ_w=(uint64_t*) circ_array;
        register uint64x1_t 	a,b,c,temp,temp1;
        register int 		i=0,count=0;
        int 			size_l=size/8,off_n;
        register int* 		ptr_circ=&off_n;

	switch(off&0x7){
		case 0: 
			off_n=off>>3;
			while(*ptr_circ<size_l){
				c=veor_u64(vld1_u64(&ref_w[i++]),vld1_u64(&circ_w[(*ptr_circ)++]));
				temp1=vadd_u64 (c,temp1);
				count=count+setbits(c);
			}
			*ptr_circ=0;
			while(i<size_l){
				c=veor_u64(vld1_u64(&ref_w[i++]),vld1_u64(&circ_w[(*ptr_circ)++]));
                                count=count+setbits(c);
			}
			
			return (size*8-count);
		case 2:
			off_n=off>>3;
			a=vld1_u64(&circ_w[(*ptr_circ)++]);
			while(*ptr_circ<size_l){
				b=vld1_u64(&circ_w[(*ptr_circ)++]);
				temp=vorr_u64(vshr_n_u64(a,16),vshl_n_u64(b,48));
				c=veor_u64(vld1_u64(&ref_w[i++]),temp);
				count=count+setbits(c);
				a=b;
			}
			*ptr_circ=0;
                        while (i<size_l){		
				b=vld1_u64(&circ_w[(*ptr_circ)++]);
                                temp=vorr_u64(vshr_n_u64(a,16),vshl_n_u64(b,48));
                                c=veor_u64(vld1_u64(&ref_w[i++]),temp);
                                count=count+setbits(c);
				a=b;
			}
			return (size*8-count);
		case 4:
			off_n=off>>3;
                        a=vld1_u64(&circ_w[(*ptr_circ)++]);
                        while(*ptr_circ<size_l){
                                b=vld1_u64(&circ_w[(*ptr_circ)++]);
                                temp=vorr_u64(vshr_n_u64(a,32),vshl_n_u64(b,32));
                                c=veor_u64(vld1_u64(&ref_w[i++]),temp);
                                count=count+setbits(c);
				a=b;
                        } 
                        *ptr_circ=0;
                        while(i<size_l){ 
                                b=vld1_u64(&circ_w[(*ptr_circ)++]);
                                temp=vorr_u64(vshr_n_u64(a,32),vshl_n_u64(b,32));
                                c=veor_u64(vld1_u64(&ref_w[i++]),temp);
                                count=count+setbits(c);
				a=b;
                        }
			return (size*8-count);
		case 6:
			off_n=off>>3;
                        a=vld1_u64(&circ_w[(*ptr_circ)++]);
                        while(*ptr_circ<size_l){
                                b=vld1_u64(&circ_w[(*ptr_circ)++]);
                                temp=vorr_u64(vshr_n_u64(a,48),vshl_n_u64(b,16));
                                c=veor_u64(vld1_u64(&ref_w[i++]),temp);
                                count=count+setbits(c);
                                a=b;    
                        }       
                        *ptr_circ=0;
                        while(i<size_l){ 
                                b=vld1_u64(&circ_w[(*ptr_circ)++]);
                                temp=vorr_u64(vshr_n_u64(a,48),vshl_n_u64(b,16));
                                c=veor_u64(vld1_u64(&ref_w[i++]),temp);
                                count=count+setbits(c);
                                a=b;
                        }
                        return (size*8-count);
		default:
			printf("something gone wrong");
			return 0;
	}	
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

