
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

	bitset<64> t(ret);

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

	Vec<GF2X> complement_arr = complement(y, m, k);

	for (int i = 0; i < complement_arr.length(); i++)
	{
		//cout << complement_arr[i] << endl;
		//cout << "a::" << a << endl;
		GF2X z = MulMod(InvMod(a, f), (complement_arr[i] - b), f);

		long dec_z = GF2X_to_Decimal(z, m);

		bool search_res = binary_search(Set.begin(), Set.end(), dec_z);

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
			//cout << "xΪ:" << x << endl;

			return 0;
		}
	}
	return 1;
}

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

			/*GF2X a = random_GF2X(m);
			GF2X b = random_GF2X(m);
			GF2X y = random_GF2X(k);*/
			GF2X a = vec_a[i];
			GF2X b = vec_b[i];
			GF2X y = vec_y[i];


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
	Hash_Table(u64 size);
	u64 Hash_fun(Elemtype value);
	void Create_HashTable(vector<long>, u64 size);
	bool unique_Lnode(Elemtype value); 
	void Insert_Lnode(Elemtype value);
	void Debug();
	vector<long> get_data(u64 index);
	u64 get_size();

private:
	Lnode* table;
	u64 num;
};

Hash_Table::Hash_Table(u64 size)
{
	this->num = size;
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
	PRNG prng(_mm_set_epi32(4253464, 3434564, 234434, 23987044));
	auto myHashSeed = prng.get<block>();
	//cout << myHashSeed << endl;
	u8 hashOut[RandomOracle::HashSize];
	RandomOracle ro_hash;
	ro_hash.Update(value);
	ro_hash.Final(hashOut);
	u64& idx = *(u64*)hashOut;
	idx %= this->num;
	return idx;
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
		else
		{
			bool reslut = this->unique_Lnode(value[i]);
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

bool Hash_Table::unique_Lnode(Elemtype value)
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
	else
	{
		bool reslut1 = this->unique_Lnode(value);
		if (reslut1)
		{
			Lnode* node = new Lnode;
			node->data = value; node->key = reslut; node->next = table[reslut].next; table[reslut].next = node;//ǰ�巨
		}
		else
			return;
	}
}

void Hash_Table::Debug()
{

	for (u64 j = 0; j < this->num; j++)
	{
		Lnode* temp = this->table[j].next;

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

void GS_Protocol_hash_search(int i, vector<long>& correctSet, Hash_Table HT, int N, int m, int k, GF2X a, GF2X b, GF2X y, GF2X f,int m2) {
	//std::cout << "-------Round<" << i << ">--------" << endl;
	//int res = 0;
	Vec<GF2X> complement_arr = complement(y, m2, k);

	for (int i = 0; i < complement_arr.length(); i++)
	{
		//cout << complement_arr[i] << endl;
		//cout << "a::" << a << endl;
		GF2X z = MulMod(InvMod(a, f), (complement_arr[i] - b), f);
		int mod = pow(2,m2);
		long dec_z = GF2X_to_Decimal(z, m) % mod;

		vector<long> res = HT.get_data(dec_z);
		for (int i = 0; i < res.size(); i++)
		{
			//cout << res[i] << endl;
			correctSet.push_back(res[i]);
		}


		//cout << "--------------------------" << endl;
	}
	//return res;
}
int main() {
	//PRNG prng(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
	//block ce = prng.get<block>();

	auto n = thread::hardware_concurrency();
	SetNumThreads(8);
	struct timespec start, finish;
	double elapsed;
	cout << n << endl;
	int Round = 5000;
	int N = pow(2, 20);
	float p = (float)1 / 2;
	//cout << "p=" << p << endl;
	float threshold = 0.5;
	int k = log2((int)(threshold * N)) + 1;
	int m = 64;
	cout << "N=" << N << endl;
	cout << "k=" << k << endl;
	cout << "m=" << m << endl;
	string f_filename = "f20.txt";
	string P1_set_filename = "p1_20.txt";
	string P1_sort_filename = "p1_sort_20.txt";
	string P2_set_filename = "p2_20.txt";
	string P2_sort_filename = "p2_sort_20.txt";
	string a_filename = "a20.txt";
	string b_filename = "b20.txt";
	string y_filename = "y20.txt";


	GF2X f = BuildIrred_GF2X(m);

	cout << "deg(f):" << deg(f) << endl;



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


	Vec<GF2X> P2_Set;
	P2_Set.SetLength(N);

	float frac = 0.1;
	int sameNum = (int)(frac * N);
	//cout << "sameNum=" << sameNum << endl;

	for (int i = 0; i < sameNum; i++) {
		P2_Set[i] = P1_Set[i];
	}

	for (int i = sameNum; i < N; i++) {

		GF2X x = random_GF2X(m);
		P2_Set.put(i, random_GF2X(m));
	}

	vector<long> P2_data;
	for (int i = 0; i < N; i++) {

		P2_data.push_back(GF2X_to_Decimal(P2_Set[i], m));
	}
	sort(P2_data.begin(), P2_data.end());

	vector<long> P1_data1;
	for (int i = 0; i < N; i++) {

		P1_data1.push_back(i);
	}

	vector<long> P2_data1;
	for (int i = 0; i < N; i++) {

		P2_data1.push_back(i);
	}


	Vec<GF2X> vec_a;
	vec_a.SetLength(Round);

	Vec<GF2X> vec_b;
	vec_b.SetLength(Round);
	Vec<GF2X> vec_y;
	vec_y.SetLength(Round);

	for (int i = 0; i < Round; i++)
	{
		vec_a.put(i, random_GF2X(m));
		vec_b.put(i, random_GF2X(m));
		vec_y.put(i, random_GF2X(k));
	}


	int m2 = k + 4;
	int size = pow(2,m2);
	Hash_Table P1_HT(size);
	P1_HT.Create_HashTable(P1_data, N);
	//P1_HT.Debug();

	Hash_Table P2_HT(size);
	P2_HT.Create_HashTable(P2_data, N);
	//P2_HT.Debug();

	int count = 0;
	int P1_count = 0;
	int P2_count = 0;
	long n1 = Round;
	//long first = 0;
	//long last = Round;
	//clock_t start = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
	//#pragma omp parallel for

	//test(n1, P1_data, P2_data, vec_a, vec_b, vec_y, N, m, k, f);

	NTL_EXEC_RANGE(n1, first, last)
		for (int i = first; i < last; i++)
		{

			/*GF2X a = random_GF2X(m);
			GF2X b = random_GF2X(m);
			GF2X y = random_GF2X(k);*/
			GF2X a = vec_a[i];
			GF2X b = vec_b[i];
			GF2X y = vec_y[i];

			vector<long> P1_correct_set, P2_correct_set;
			GS_Protocol_hash_search(i, P1_correct_set, P1_HT, N, m, k, a, b, y, f,m2);
			GS_Protocol_hash_search(i, P2_correct_set, P2_HT, N, m, k, a, b, y, f,m2);

			if (compare2(P1_correct_set, P2_correct_set))
			{
				count += 1;
			}
		}
	NTL_EXEC_RANGE_END

		float ex = (float)5 / 8;
		int target = (int)(ex * p * Round);
		cout << "count =" << count << endl;
		cout << "target =" << target << endl;
		//clock_t end = clock();


		clock_gettime(CLOCK_MONOTONIC, &finish);

	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	cout << "Total time===" << elapsed << endl;
	return 0;
}