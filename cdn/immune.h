#include "config.h"
#include "define.h"
#include <iostream>
#include <vector>
#include <memory.h>
#include <algorithm>
using namespace std;
//#define _GENERATION_AMOUNT 20
//定义全局常量
const int C = _PONISH_FACTOR; //ponishment
const float ALPHA = _ALPHA; //im_Probability=ALPHA*im_affinity+(1-ALPHA)*im_density

typedef vector<uint_16> T;
typedef vector<float> F;
//typedef vector<int> P;
//typedef vector<float> F;
typedef struct _Map_im_Probability		//排序映射表
{
	int i;
	float vProbability;
}Map_imProbability;

bool comp1(const Map_imProbability &a,const Map_imProbability &b)	//排序用函数
{
    return a.vProbability>b.vProbability;
}

bool comp2(const float &a,const float &b)	//排序用函数
{
    return a>b;
}

inline int fnRndBoundary(int iBegin, int iEnd)	//产生随机的下标
{

	return rand()%(iEnd-iBegin) + iBegin;
}




typedef struct Immune_global_t
{
  int im_number;        //表示未知数个数
  int im_maxserver;
  int im_scost;
  int deno;             //表示im_number*_GENERATION_AMOUNT
  float pcross;
  float pmutation;
}Immune_global;
//需要定义T,P,STAT
//T是群体基因，P是适应度, F表示浮点向量
class immune
{
public:
  immune();
  immune(FLOW flow, CUSTOMER customer, uint_16 cost);    //OK构造函数（两个参数都是vector数组）
  ~immune();

  //bool imGetCost(const int cost); //获取每台服务器成本
  vector<T> imSendSpecies();       //把产生的最新种群发送给最小费用最大流
  bool imReceiveAfinity();    //获取最小费用最大流产生的亲和度信息,包括链路费用 和 可行性

  bool imCreateGene();        //产生初始种群
  bool imGeneChoose();        //基因选择
  bool imGeneMix();           //基因交叉
  bool imGeneAberrance();     //基因变异
  bool imGeneSort();			//期望繁殖率排序
  bool im_memoryGenneSort();	//更新记忆库
	bool im_RETGen_Gene();		//返回新的种群

  //bool imEvalAll();           //测试所有基因适应度
  bool imAfinity(Affinity_Info cost_err);           //测试所有基因亲和度
  //bool imAfinityOne(T &Gene);    //测试某一个基因亲和度
  bool imDensity();           //测试种群密度/浓度
  int imDensityOne(T &Gene);    //测试某一个基因同其他基因的相似度
  bool imBreedProbability();    //繁殖概率

  //bool imPutGenetoMemory(T Gene);   //把一个基因放入记忆细胞中
  //bool imGetGenefromMemory();   //从记忆细胞中获取基因

  void imPrintInfo();

  vector<T> im_Gen_Gene;  //种群信息
  vector<T> im_memory;  //记忆细胞
  vector<T> im_Gen_Gene_Father;   //父代种群

  F im_Probability;
  F im_memoryProbability;
  F im_affinity;
  F im_density;
  //float im_Probability[_GENERATION_AMOUNT];     //P基因选择概率
  //float im_memory_Probability[_GENERATION_AMOUNT];  //记忆库基因选择概率
  //float im_affinity[_GENERATION_AMOUNT];        //F亲和度
  //float im_density[_GENERATION_AMOUNT];         //F浓度

  Immune_global ginfo;  //全局信息
private:
  bool imGeneAberranceOne(int index);   //变异某个基因
  //bool imWheel(T Gene);     //轮盘赌
};

immune::immune()
{
}


immune::~immune()
{
}

immune::immune(FLOW flow, CUSTOMER customer, uint_16 cost)
{//初始化全局变量
  //返回可行解放在第一行
  int c_num = customer.size();
  int n_num = flow.size()-c_num;
  ginfo.im_number = n_num;
  ginfo.im_scost = cost;
  ginfo.im_maxserver = c_num;
  ginfo.deno = ginfo.im_number * _GENERATION_AMOUNT;   //可以放在构造函数中优化
  ginfo.pcross = _PCROSS;		//交叉概率
  ginfo.pmutation = _PMUTATION;	//变异概率
  //将im_Gen_Gene和im_memory初始化为0

  im_Gen_Gene.reserve(_GENERATION_AMOUNT);  //种群信息
  im_memory.reserve(_GENERATION_AMOUNT);  //记忆细胞
  im_Gen_Gene_Father.reserve(_GENERATION_AMOUNT);   //父代种群
  im_Probability.reserve(_GENERATION_AMOUNT);
  im_memoryProbability.reserve(_GENERATION_AMOUNT);
  im_affinity.reserve(_GENERATION_AMOUNT);
  im_density.reserve(_GENERATION_AMOUNT);
  int i, j;
  T temp(ginfo.im_number, 0);
  for(i = 0;i < _GENERATION_AMOUNT;i++)
  {
    im_Gen_Gene.push_back(temp);
    im_Probability.push_back(0);
    im_affinity.push_back(0);
    im_density.push_back(0);
  }

  //im_Gen_Gene[i].assign(flsize, 0);
  for(j = 0;j < c_num;j++)
  {
    for(i = 0;i < n_num;i++)
    {
      if(flow[i][j+n_num] > 0)
        im_Gen_Gene[0][i] = 1;
        continue;
    }
  }
}

bool immune::imCreateGene()        //随机产生初始种群
{
	int i,j,temp,cnt=0,rand_count,size;
	srand((unsigned int)time(0));
	T first=im_Gen_Gene[0];
	T *second;
	size=first.size();
	T::iterator it;
  cnt = count(first.begin(), first.end(), 1);
	/*for(it=first.begin();it!=first.end();it++)
	{
		if((*it)==1)
			count++;
	}*/
	int MIN_SOLUTION = cnt/2;
	for(i=1;i<_GENERATION_AMOUNT;i++)
	{
		second=&im_Gen_Gene[i];
		//im_Gen_Gene[i].assign(size,0);					//先将该种群个体初始化清零
		rand_count=rand()%(cnt-MIN_SOLUTION)+MIN_SOLUTION;
		for(j=0;j<rand_count;j++)
		{
			temp=rand()%size;
			if((*second)[temp]==1)
			{
				j--;
				continue;
			}
			else
				(*second)[temp]=1;
		}
	}
  return true;
}

bool immune::imGeneChoose()        //基因选择
{
		int i,j,size;
		imGeneSort();		//先对当前的种群进行一个降序排序，排序依据为期望繁殖概率

		for(i=0;i<im_N;i++)	//取前im_N个作为父代群体
		{
			im_Gen_Gene_Father.push_back(im_Gen_Gene[i]);
		}

		if(im_memory.empty())
		{
			for(i=0;i<im_m;i++)	//取前im_m个放入到记忆库中
			{
				im_memory.push_back(im_Gen_Gene[i]);
				im_memoryProbability.push_back(im_Probability[i]);
			}
		}
		else					//记忆库非空，此时需要放入后更新记忆库
		{
			size = im_memory.size();
			if((_GENERATION_AMOUNT-size)<im_m)
			{
				for(i = 0,j= _GENERATION_AMOUNT-im_m-1;i<im_m && j<_GENERATION_AMOUNT; i++,j++)
				{
					im_memory[j] = im_Gen_Gene[i];
					im_memoryProbability[j] = im_Probability[i];
				}
			}
			else
			{
				for(i = 0;i<im_m; i++)
				{
					im_memory.push_back(im_Gen_Gene[i]);
					im_memoryProbability.push_back(im_Probability[i]);
				}
			}

			im_memoryGeneSort();		//更新记忆库
			sort(im_memoryProbability.begin(),im_memoryProbability.end(),comp2);	//更新记忆库对应的概率
		}
    return true;
}

bool immune::imGeneSort()								//种群排序，按照期望繁殖概率
{
	int i;
	Map_imProbability *p = new Map_imProbability[_GENERATION_AMOUNT];
	for(i=0;i<_GENERATION_AMOUNT;i++)
	{
		p[i].i = i;
		p[i].vProbability = im_Probability[i];
	}
	sort(p,p+_GENERATION_AMOUNT,comp1);
	vector<T> temp;			//辅助vector
	for(i=0;i<_GENERATION_AMOUNT;i++)
	{
		temp[i]=im_Gen_Gene[p[i].i];
	}
	for(i=0;i<_GENERATION_AMOUNT;i++)	//将排序后的种群复制到初始种群中
	{
		im_Gen_Gene[i] = temp[i];
	}
	delete p;
	p = NULL;
	sort(im_Probability.begin(),im_Probability.end(),comp2);	//对相应的种群期望概率进行排序
  return true;
}

bool immune::im_memoryGeneSort()
{
	int i,j,size;
	Map_imProbability *p = new Map_imProbability[_GENERATION_AMOUNT];
	for(i=0;i<_GENERATION_AMOUNT;i++)
	{
		p[i].i = i;
		p[i].vProbability = im_memoryProbability[i];
	}
	sort(p,p+_GENERATION_AMOUNT,comp1);
	vector<T> temp;			//辅助vector
	for(i=0;i<_GENERATION_AMOUNT;i++)
	{
		//temp[i].assign(im_memory[p[i].i].begin(),im_memory[p[i].i].end());
		temp[i] = im_memory[p[i].i];
	}
	for(i=0;i<_GENERATION_AMOUNT;i++)	//将排序后的记忆基因更新至记忆库中
	{
		//im_memory[i].assign(temp[i].begin(),temp[i].end());
		im_memory[i] = temp[i];
	}
	delete p;
	p = NULL;
}

bool immune::imGeneMix()
{
	srand((unsigned int )time(0));
	//std::vector<float> temp=im_affinity;	//temp：临时保存亲和度
  vector<T> temp(_GENERATION_AMOUNT, T(ginfo.im_number, 0));			//辅助的父代种群
	double pick;							//随机产生的交叉概率
	int i,size;
	int iFather;                         //父亲的代号
	int iMother;                         //母亲的代号
	T FatherBK(_GENERATION_AMOUNT);
  T MotherBK(_GENERATION_AMOUNT);				//父母的基因
	T Child1(_GENERATION_AMOUNT);
  T Child2(_GENERATION_AMOUNT);                     //父亲与母亲杂交出的子女的基因
	T::iterator V_iter;
	int Low_index;
	int High_index;
	int im_Dvalue;
	size = im_Gen_Gene_Father[0].size();
	for(i=0; i<_P_GENE_MIX; i++)
	{
		pick = rand()/(double)(RAND_MAX);
		if(pick>ginfo.pcross)
			continue;
		iFather=_ITEMP;
		do
		{
			iMother = _ITEMP;
		}while(iMother == iFather);
		//Child1.reserve(size);         //初始化子女的碱基数
		//Child2.reserve(size);
		//Child1.clear();
		//Child2.clear();
		Low_index = fnRndBoundary(0, im_N-2);
		High_index = fnRndBoundary(Low_index+1, im_N-1);

		FatherBK = im_Gen_Gene_Father[iFather];
		MotherBK = im_Gen_Gene_Father[iMother];
		copy (FatherBK.begin()+Low_index, FatherBK.begin()+High_index+1,\
				   back_inserter(Child1));

		copy (MotherBK.begin()+Low_index, MotherBK.begin()+High_index+1,\
				   back_inserter(Child2));

		rotate (FatherBK.begin(), FatherBK.begin()+High_index+1, FatherBK.end());
		rotate (MotherBK.begin(), MotherBK.begin()+High_index+1, MotherBK.end());

		for (V_iter = im_Gen_Gene_Father[iFather].begin()+Low_index;\
				V_iter != im_Gen_Gene_Father[iFather].begin()+High_index+1; ++V_iter)
		{
			MotherBK.erase(std::remove(MotherBK.begin(), MotherBK.end(), *V_iter),\
						   MotherBK.end());
		}

		for (V_iter = im_Gen_Gene_Father[iMother].begin()+Low_index;\
				V_iter != im_Gen_Gene_Father[iMother].begin()+High_index+1; ++V_iter)
		{

			FatherBK.erase(std::remove(FatherBK.begin(), FatherBK.end(), *V_iter),\
						   FatherBK.end());
		}
		im_Dvalue = size -High_index - 1;
		copy(MotherBK.begin(), MotherBK.begin()+im_Dvalue, back_inserter(Child1));
		copy(MotherBK.begin()+im_Dvalue, MotherBK.end(), inserter(Child1,Child1.begin()));

		copy(FatherBK.begin(), FatherBK.begin()+im_Dvalue, back_inserter(Child2));
		copy(FatherBK.begin()+im_Dvalue, FatherBK.end(), inserter(Child2,Child2.begin()));

		temp[2*i]  = Child1;
		temp[2*i+1] = Child2;
	}
	for(i=0; i<_GENERATION_AMOUNT; i++)
	{
		im_Gen_Gene_Father[i] = temp[i];
	}
  return true;
}

bool immune::imGeneAberrance()     //基因变异
{
	srand((unsigned int)time(0));
	int i;
	double RVariation;						//随机的变异概率
	for(i=0; i<im_N; i++)
	{
		RVariation = rand()/(double)(RAND_MAX);
		if(RVariation > ginfo.pmutation)		//发生变异
		{
			imGeneAberranceOne(i);
		}
	}
	im_RETGen_Gene();					//更新新一代种群
  return true;
}

bool immune::im_RETGen_Gene()  //返回新的种群，根据交叉变异后的父代和记忆库中的前m个记忆基因
{
	int i,j;
	for(i=0,j=0; i<_GENERATION_AMOUNT; i++)
	{
		if(i<im_N)
		{
			im_Gen_Gene[i] = im_Gen_Gene_Father[i];
		}
		else
		{
			im_Gen_Gene[i] = im_memory[j];
			j++;
		}
	}
  return true;
}

bool immune::imGeneAberranceOne(int index)  //变异某个基因
{
	srand((unsigned int)time(0));
	int temp;
	int size = im_Gen_Gene_Father[0].size();
	T::iterator V_it;
	//int Low_limit = 0;
	//int High_limit = size-1;
	temp = rand()%size;
	while(temp == size)
		temp = rand()%size;
	int num = count(im_Gen_Gene_Father[index].begin(),im_Gen_Gene_Father[index].end(),1);
	if(im_Gen_Gene_Father[index][temp] == 1)
		im_Gen_Gene_Father[index][temp] = 0;
	else
	{
		if(num>=ginfo.im_maxserver)
		{
			while(im_Gen_Gene_Father[index][temp] == 0)
			{
				temp = rand()%size;
				while(temp == size)
					temp = rand()%size;
			}
			im_Gen_Gene_Father[index][temp] = 0;
		}
		else
			im_Gen_Gene_Father[index][temp] = 1;
	}
}

bool immune::imAfinity(Affinity_Info cost_err)          //测试所有基因亲和度
{//亲和度由最小费用最大流传入的 链路总成本 和 是否是可行解 计算得出
  int i;
  float probability;
  //测试每组基因的亲和度
  for(i = 0;i < _GENERATION_AMOUNT; i++)
  {
    probability = count(im_Gen_Gene[i].begin(),im_Gen_Gene[i].end(),1)*ginfo.im_scost\
                + cost_err[i].link_cost - C*cost_err[i].flow_err; //定义C
    probability = (probability>0) ? (1/probability):0;
    im_affinity[i] = probability;
  }
  return true;
}

bool immune::imDensity()           //测试种群密度/浓度
{
  int i;
  for(i = 0;i < _GENERATION_AMOUNT; i++)
  {
    im_density[i] = (ginfo.deno-imDensityOne(im_Gen_Gene[i]))/ginfo.deno;
  }
  return true;
}

int immune::imDensityOne(T &Gene)    //测试某一个基因同其他基因的相似度
{//通过异或操作判断基因之间的相似度
  int i, j;
  int count = 0;
  for(i = 0;i < _GENERATION_AMOUNT; i++)
    for(j = 0;j < ginfo.im_number;j++)
      count += (Gene[j]^im_Gen_Gene[i][j]);
  return count;
}

bool immune::imBreedProbability()
{
  int i;
  float sum_afinity=0, sum_dencity=0;
  for(i = 0;i < _GENERATION_AMOUNT;i++)
  {
    sum_afinity += im_affinity[i];
    sum_dencity += im_density[i];
  }
  float aff, den;
  for(i = 0;i < _GENERATION_AMOUNT;i++)
  {
    aff = im_affinity[i]/sum_afinity;
    den = im_density[i]/sum_dencity;
    im_Probability[i] = ALPHA*(aff-den)+den;    //define中定义ALPHA
  }
  return true;
}

void immune::imPrintInfo()
{//打印出
  int i, j;
  cout<<"Output Species Info"<<endl;
  for(i = 0;i < _GENERATION_AMOUNT;i++)
  {
    for(j = 0;j < ginfo.im_number;j++)
      cout<<im_Gen_Gene[i][j]<<" ";
    cout<<endl;
  }
  cout<<"Output Affinity Info"<<endl;
  for(i = 0;i < ginfo.im_number;i++)
    cout<<im_affinity[i]<<" ";
  cout<<endl;
  cout<<"Output Density Info"<<endl;
  for(i = 0;i < ginfo.im_number;i++)
    cout<<im_density[i]<<" ";
  cout<<endl;
  cout<<"Output Probability Info"<<endl;
  for(i = 0;i < ginfo.im_number;i++)
    cout<<im_Probability[i]<<" ";
  cout<<endl;
}
