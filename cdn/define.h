#ifndef __DEFINE_H__
#define __DEFINE_H__


#define _GENERATION_AMOUNT 28
#define _PCROSS 0.5
#define _PMUTATION 0.12
#define _ITEMP rand()%(((_GENERATION_AMOUNT+1)*_GENERATION_AMOUNT)/2)	//TODO：随机的下标

#define im_m (_GENERATION_AMOUNT/4)				//TODO:这样每次都要计算不合理；每次放入记忆库的个体数量
#define im_N (_GENERATION_AMOUNT/2+_GENERATION_AMOUNT/4)				//TODO：每一代父代群体的个数
#define	_P_GENE_MIX (_GENERATION_AMOUNT-1)/2  		//TODO：杂交次数
#define _PONISH_FACTOR 10
#define _ALPHA 0.7
#endif
