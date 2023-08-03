
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include "cryptoTools/Crypto/PRNG.h"
#include "cryptoTools/Crypto/Commit.h"
#include "cryptoTools/Common/Log.h"
#include <NTL/GF2X.h>
#include <NTL/GF2.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>
#include <omp.h>
#include<thread>
#include<vector>
#include<ctime>
#include<cmath>
#include <sstream>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <time.h>
#include "sha1.hpp"
using namespace boost;
using namespace std;
using namespace NTL;
using namespace osuCrypto;

typedef u64 Elemtype;

string GF2X_to_String(GF2X x) {
	vector<char> arr;
	for (int i = 0; i < deg(x) + 1; i++)
	{
		if (x[i] == GF2(0))
		{
			arr.push_back('0');
		}
		else
		{
			arr.push_back('1');
		}
	}
	string ret;
	for (int i = 0; i < arr.size(); i++)
	{
		ret += arr[i];
	}
	//cout << "转换后的bit串为:" << ret << endl;
	//bitset<64> t(ret);
	//cout << "转换为整型为:" << t.to_ulong() << endl;
	return ret;
}

unsigned long GF2X_to_Decimal(GF2X x, int m) {
	vector<char> arr;
	for (int i = 0; i < deg(x) + 1; i++)
	{
		if (x[i] == GF2(0))
		{
			arr.push_back('0');
		}
		else
		{
			arr.push_back('1');
		}
	}
	string ret;
	for (int i = 0; i < arr.size(); i++)
	{
		ret += arr[i];
	}
	for (int i = 0; i < m - deg(x) - 1; i++)
	{
		ret += '0';
	}
	//cout << "原来的x为:" << x << endl;
	//cout << "转换后的bit串为:" << ret << endl;
	//dynamic_bitset<> db(m, ret);
	//cout << "转换后的bit串为:" << ret << endl;
	bitset<64> t(ret);
	//cout << "转换为整型为:" << t.to_ulong() << endl;
	return t.to_ulong();
}

vector<string> generateSequence(const int n) {

	vector<string> vec_Sequence;
	//bitset b = new bitset(n);


	for (int i = 0; i < pow(2, n); i++)
	{
		dynamic_bitset<> db(n, i);
		string s;
		to_string(db, s);
		//cout << bitset<n>(i) << endl;
		vec_Sequence.push_back(s);
	}
	return vec_Sequence;
}



Vec<GF2X> complement(GF2X x, int m, int k) {
	Vec<GF2X> complement_arr;
	//int n = m - 1 - deg(x);
	int n = m - k;
	vector<string> vec_Sequence = generateSequence(n);

	for (int i = 0; i < vec_Sequence.size(); i++)
	{
		GF2X a = x;
		for (int j = k, l = 0; j < m - 1, l < vec_Sequence[i].length(); j++, l++)
		{
			SetCoeff(a, j, vec_Sequence[i][l]);
		}
		complement_arr.append(a);
	}
	return complement_arr;
}

GF2X GS_Hash(GF2X a, GF2X b, GF2X x, int m, int k, GF2X f) {
	//GF2X f = GF2X();
	//f.SetLength(k);
	//f.normalize();
	//GF2X hash_result = (MulTrunc(a, x, m) + b);

	GF2X result = MulMod(a, x, f) + b;
	GF2X hash_result;

	for (int i = 0; i < k; i++)
	{
		SetCoeff(hash_result, i, result[i]);
	}
	//cout << "result" << result << endl;
	//cout << "hash_result>" << hash_result << endl;
	//cout << "deg(hash_result)" << deg(hash_result) << endl;
	return hash_result;
}



int GS_Protocol_sort_search(int i, vector<long>& correctSet, vector<long> Set, int N, int m, int k, GF2X a, GF2X b, GF2X y, GF2X f) {
	//std::cout << "-------Round<" << i << ">--------" << endl;
	int res = 0;
	//对y进行补位
	Vec<GF2X> complement_arr = complement(y, m, k);
	// 二分查找判断
	for (int i = 0; i < complement_arr.length(); i++)
	{
		//cout << complement_arr[i] << endl;
		//cout << "a::" << a << endl;
		GF2X z = MulMod(InvMod(a, f), (complement_arr[i] - b), f);
		//cout << "当前y=" << y << endl;
		//cout << "补位y=" << complement_arr[i] << endl;
		//cout << "当前z=" << z << endl;
		long dec_z = GF2X_to_Decimal(z, m);
		//cout << "转化z=" << dec_z << endl;
		bool search_res = binary_search(Set.begin(), Set.end(), dec_z);
		//cout << "是否找到=" << search_res << endl;
		if (search_res)
		{
			//cout << dec_z << endl;
			correctSet.push_back(dec_z);
			res = 1;
		}
		//cout << "--------------------------" << endl;
	}
	return res;
}
int GS_Protocol(int i, Vec<GF2X>& correctSet, Vec<GF2X> Set, int N, int m, int k, GF2X a, GF2X b, GF2X y, GF2X f) {
	//std::cout << "-------Round<" << i << ">--------" << endl;
	int res = 0;
	//#pragma omp parallel for
	for (int i = 0; i < N; i++) {

		GF2X x = Set[i];
		GF2X hash_value = GS_Hash(a, b, x, m, k, f);
		if (hash_value == y)
		{
			//std::cout << i << endl;
			correctSet.append(x);
			res = 1;
			// 测试
			//cout << "当前的a值为:" << a << endl;
			//cout << "当前的b值为:" << b << endl;
			//cout << "当前的y值为:" << y << endl;

			// //对y进行补位
			//Vec<GF2X> complement_arr = complement(y, m,k);
			//for (int i = 0; i < complement_arr.length(); i++)
			//{
			//	//cout << "补位:" << complement_arr[i] << endl;
			//	GF2X z = MulMod(InvMod(a, f), (complement_arr[i] - b), f);
			//	//GF2X h = GS_Hash(a, b, z, m, k, f);
			//	//cout << "当前的h值为:" << h << endl;
			//	//cout << "计算出的z:" << z << endl;
			//	//cout << "当前的的x:" << x << endl;
			//	if (z == x)
			//	{
			//		cout << "********找到********" << endl;
			//	}
			//}

		}
	}
	return res;
}
int compare(Vec<GF2X> setA, Vec<GF2X> setB) {
	for (int i = 0; i < setA.length(); i++)
	{
		for (int j = 0; j < setB.length(); j++)
		{
			if (setA[i] == setB[j])
			{
				return 1;
			}
		}
	}
	return 0;
}
int compare2(vector<long> setA, vector<long> setB) {
	for (int i = 0; i < setA.size(); i++)
	{
		for (int j = 0; j < setB.size(); j++)
		{
			if (setA[i] == setB[j])
			{
				return 1;
			}
		}
	}
	return 0;
}
int judge(Vec<GF2X> set, GF2X x) {

	if (set.position(x) == -1)
	{
		return 1;
	}
	return 0;
}

int judge2(Vec<GF2X> set, GF2X x) {

	for (int i = 0; i < set.length(); i++)
	{
		//cout <<"set:::" <<set[i] << endl;
		//cout << "x:::" << x << endl;
		//bool res = set[i] == x;
		//cout << "res:::" <<res  << endl;
		if (set[i] == x)
		{
			//cout << "x为:" << x << endl;
			//cout << "相等的set为:" << set[i] << endl;
			return 0;
		}
	}
	return 1;
}
//模板函数：将string类型变量转换为常用的数值类型（此方法具有普遍适用性） 
template <class Type>
Type stringToNum(const string& str) {
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}


void test(long n1, vector<long> P1_data, vector<long> P2_data, Vec<GF2X> vec_a, Vec<GF2X> vec_b, Vec<GF2X> vec_y, int N, int m, int k, GF2X f) {
	int count = 0;
	int P1_count = 0;
	int P2_count = 0;
	NTL_EXEC_RANGE(n1, first, last)
		for (int i = first; i < last; i++)
		{
			// 从GF(2^n)中随机选择a、b、y
			/*GF2X a = random_GF2X(m);
			GF2X b = random_GF2X(m);
			GF2X y = random_GF2X(k);*/
			GF2X a = vec_a[i];
			GF2X b = vec_b[i];
			GF2X y = vec_y[i];


			//// 初始
			//Vec<GF2X> P1_correct_set, P2_correct_set;
			//P1_count += GS_Protocol(i, P1_correct_set, P1_Set, N, m, k, a, b, y, f);
			//P2_count += GS_Protocol(i, P2_correct_set, P2_Set, N, m, k, a, b, y, f);


			// 优化
			vector<long> P1_correct_set, P2_correct_set;
			P1_count += GS_Protocol_sort_search(i, P1_correct_set, P1_data, N, m, k, a, b, y, f);
			P2_count += GS_Protocol_sort_search(i, P2_correct_set, P2_data, N, m, k, a, b, y, f);


			if (compare2(P1_correct_set, P2_correct_set))
			{
				count += 1;
			}
		}
	NTL_EXEC_RANGE_END

}



class Lnode
{
public:
	bool IsEmpty()
	{
		if (this->next == nullptr)
			return true;
		else
			return false;
	}
public:
	Lnode* next;
	u64 key;
	Elemtype data;
};


class Hash_Table
{

public:
	Hash_Table(u64 size,block seed);
	u64 Hash_fun(Elemtype value);//定义哈希函数
	void Create_HashTable(vector<long>, u64 size);//创建哈希表
	bool unique_Lnode(Elemtype value); //判断结点是否重复
	void Insert_Lnode(Elemtype value);//插入结点进入哈希表
	void Debug();//测试一下
	~Hash_Table();
	vector<long> get_data(u64 index);
	u64 get_size();

private:
	Lnode* table;
	u64 num;
	block seed;
};

Hash_Table::Hash_Table(u64 size,block seed)
{
	this->num = size;
	this->seed = seed;
	table = new Lnode[num];
	for (u64 i = 0; i < this->num; i++)
	{
		this->table[i].data = NULL;
		this->table[i].key = i;
		this->table[i].next = nullptr;
	}
}

u64 Hash_Table::Hash_fun(Elemtype value)
{

	//cout << this->seed << endl;
	PRNG prng(this->seed);
	auto myHashSeed = prng.get<block>();
	RandomOracle ro_hash;
	//cout << myHashSeed << endl;
	vector<RandomOracle> mHashs;
	mHashs.resize(1);
	mHashs[0].Update(prng.get<block>());
	auto hash1 = mHashs[0];
	//auto hash1 = ro_hash.Update(prng.get<block>());
	u8 hashOut[RandomOracle::HashSize];
	
	hash1.Update(value);
	hash1.Final(hashOut);
	u64& idx = *(u64*)hashOut;
	idx %= this->num;
	return idx;
}
Hash_Table::~Hash_Table()//尤其要注意开辟在堆区的数据的删除
{
	for (u64 j = 0; j < this->num; j++)
	{
		Lnode* p1, * p2;
		p1 = this->table[j].next;
		p2 = p1;
		while (p1)
		{
			p2 = p1->next;
			delete p1;
			p1 = p2;
		}
		this->table[j].next = nullptr;
	}
	delete[] this->table;
	//cout << "释放空间成功！" << endl;
}
void Hash_Table::Create_HashTable(vector<long> value, u64 size)
{
	for (u64 i = 0; i < size; i++)
	{
		//cout << i << endl;

		u64 res = this->Hash_fun(value[i]);
		//cout << res << endl;
		//cout << "---" << endl;
		u64 _key = this->Hash_fun(value[i]);

		if (table[_key].IsEmpty())
		{
			Lnode* node = new Lnode;
			node->data = value[i]; node->key = _key; node->next = nullptr; table[_key].next = node;
		}
		else//如果现在结点不为空，先去遍历，看是否该元素已经添加
		{
			bool reslut = this->unique_Lnode(value[i]);//reslut为1：该元素为新元素，为0，就不再添加了
			if (reslut)
			{
				Lnode* node = new Lnode;
				node->data = value[i]; node->key = _key; node->next = table[_key].next; table[_key].next = node;
			}
			else
				continue;
		}
	}
}

bool Hash_Table::unique_Lnode(Elemtype value)//判断数据是否已经存在
{
	u64 index = this->Hash_fun(value);
	Lnode* temp = table[index].next;
	while (temp)
	{
		if (temp->data == value)
			return false;
		temp = temp->next;
	}
	return true;
}

void Hash_Table::Insert_Lnode(Elemtype value)
{
	u64 reslut = this->Hash_fun(value);
	if (this->table[reslut].next == nullptr)
	{
		Lnode* node = new Lnode;
		node->data = value; node->key = reslut; node->next = nullptr; table[reslut].next = node;
	}
	else//如果现在结点不为空，先去遍历，看是否该元素已经添加
	{
		bool reslut1 = this->unique_Lnode(value);//reslut为1：该元素为新元素，为0，就不再添加了
		if (reslut1)
		{
			Lnode* node = new Lnode;
			node->data = value; node->key = reslut; node->next = table[reslut].next; table[reslut].next = node;//前插法
		}
		else
			return;
	}
}

void Hash_Table::Debug()
{
	cout << "测试如下：" << endl;
	for (u64 j = 0; j < this->num; j++)
	{
		Lnode* temp = this->table[j].next;
		cout << "索引为" << j << "的数值：";
		if (temp == nullptr)
		{
			cout << "NULL" << endl;
			continue;
		}
		while (temp)
		{
			cout << temp->data << "\t";
			temp = temp->next;
			if (!temp)
				cout << endl;
		}
	}
}

u64 Hash_Table::get_size()
{
	u64 count = 0;
	for (u64 i = 0; i < this->num; i++)
	{
		if (this->table[i].next != nullptr)
		{
			count++;
		}
	}
	return count;
}
vector<long> Hash_Table::get_data(u64 index) {
	vector<long> result;
	Lnode* temp = this->table[index].next;
	//cout << "索引为" << j << "的数值：";
	if (temp == nullptr)
	{
		//cout << "NULL" << endl;
		return result;
	}
	while (temp)
	{
		//cout << temp->data << "\t";
		result.push_back(temp->data);
		temp = temp->next;
		//if (!temp)
		//	cout << endl;
	}
	return result;
}

void GS_Protocol_hash_search(int i, vector<long>& correctSet, Hash_Table HT, int N, int m, int k, GF2X a, GF2X b, GF2X y, GF2X f, int m2) {
	//std::cout << "-------Round<" << i << ">--------" << endl;
	//int res = 0;
	//对y进行补位
	Vec<GF2X> complement_arr = complement(y, m2, k);

	for (int i = 0; i < complement_arr.length(); i++)
	{
		//cout << complement_arr[i] << endl;
		//cout << "a::" << a << endl;
		GF2X z = MulMod(InvMod(a, f), (complement_arr[i] - b), f);
		//cout << "当前y=" << y << endl;
		//cout << "补位y=" << complement_arr[i] << endl;
		//cout << "当前z=" << z << endl;
		int mod = pow(2, m2);
		long dec_z = GF2X_to_Decimal(z, m) % mod;
		//cout << "转化后的dec_z=" << dec_z << endl;
		vector<long> res = HT.get_data(dec_z);
		for (int i = 0; i < res.size(); i++)
		{
			//cout << res[i] << endl;
			correctSet.push_back(res[i]);
		}
		//cout << "转化z=" << dec_z << endl;
		//bool search_res = binary_search(Set.begin(), Set.end(), dec_z);
		//cout << "是否找到=" << search_res << endl;

		//cout << "--------------------------" << endl;
	}
	//return res;
}
int main() {
	//PRNG prng(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
	//block ce = prng.get<block>();

	//cout << "随机种子=" << r << endl;
	
	//PRNG prng1(aa);
	//auto myHashSeed = prng1.get<block>();
	//cout << myHashSeed << endl;
	auto n = thread::hardware_concurrency();//获取cpu核心个数 
	SetNumThreads(16);
	struct timespec start, finish;
	double elapsed;
	cout << n << endl;
	int Round = 100;
	int N = pow(2, 18);
	float p = (float)1 / 2;
	//cout << "p=" << p << endl;
	float threshold = 0.5;
	int k = log2((int)(threshold * N)) + 1;
	int m = 64;
	cout << "N=" << N << endl;
	cout << "k为:" << k << endl;
	cout << "m为:" << m << endl;


	//// 写文件
	GF2X f = BuildIrred_GF2X(m);

	cout << "deg(f):" << deg(f) << endl;


	// 生成参与方P1的集合
	Vec<GF2X> P1_Set;
	P1_Set.SetLength(N);
	for (int i = 0; i < N; i++) {

		GF2X x = random_GF2X(m);
		P1_Set.put(i, random_GF2X(m));

	}

	vector<long> P1_data;
	for (int i = 0; i < N; i++) {

		P1_data.push_back(GF2X_to_Decimal(P1_Set[i], m));
	}
	sort(P1_data.begin(), P1_data.end());

	int m2 = k + 8;
	//int size = pow(2, m2);
	int size = 8 * N;
	int count = 10000;
	int flag = 0;
	vector<int> size_arr;
#pragma omp parallel for
	for (int i = 0; i < count; i++)
	{
		// 根据时间戳生成随机种子
		srand((unsigned)time(NULL));
		u64 r = rand();
		block seed = toBlock(r+i);

		// 创建hashTable
		//cout << seed << endl;
		Hash_Table P1_HT(size,seed);
		P1_HT.Create_HashTable(P1_data, N);
		//P1_HT.Debug();
		flag += P1_HT.get_size();
		size_arr.push_back(P1_HT.get_size());
		//cout << "大小:" << P1_HT.get_size() << endl;

	}

	sort(size_arr.begin(), size_arr.end());
	//for (int i = 0; i < size_arr.size(); i++)
	//{
	//	cout  << size_arr[i] << endl;
	//}
	
	for (int i = 0; i < 10; i++)
	{
		 cout << size_arr[i] << endl;
	}

	int result = flag / count;
	cout << "集合大小为:::" << N << endl;
	cout << "哈希表平均大小:::" << result << endl;


	
	return 0;
}