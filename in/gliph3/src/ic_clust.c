#include "ic_clust.h"

// const int NUM_AA_SCORE_MATRIX = 24;
const char SCORE_MATRIX_AAS[] ="ARNDCQEGHILKMFPSTWYVBZX*";
const char BLOSUM62[]="\
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\
";

const int BLOSUM62_SCORE_MATRIX[NUM_AA_SCORE_MATRIX][NUM_AA_SCORE_MATRIX]={
	{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
	{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
	{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
	{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
	{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
	{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
	{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
	{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
	{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
	{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
	{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
	{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
	{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
	{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
	{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
	{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
	{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
	{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
	{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
	{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
	{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
	{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
	{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
	{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};

int* convert_cseq_to_iseq(const char* cseq)
{
	int* iseq = icalloc(strlen(cseq)+1, int);
	int idx = -1;
	for(int i = 0; i < strlen(cseq); i++)
	{
		switch (cseq[i])
		{
		case 'A': idx = 0; break;
		case 'R': idx = 1; break;
		case 'N': idx = 2; break;
		case 'D': idx = 3; break;
		case 'C': idx = 4; break;
		case 'Q': idx = 5; break;
		case 'E': idx = 6; break;
		case 'G': idx = 7; break;
		case 'H': idx = 8; break;
		case 'I': idx = 9; break;
		case 'L': idx = 10; break;
		case 'K': idx = 11; break;
		case 'M': idx = 12; break;
		case 'F': idx = 13; break;
		case 'P': idx = 14; break;
		case 'S': idx = 15; break;
		case 'T': idx = 16; break;
		case 'W': idx = 17; break;
		case 'Y': idx = 18; break;
		case 'V': idx = 19; break;
		case 'B': idx = 20; break;
		case 'Z': idx = 21; break;
		case 'X': idx = 22; break;
		case '*': idx = 23; break;
		default: idx = 23;
		}
		iseq[i] = idx;
	}
	iseq[strlen(cseq)] = -1;
	return iseq;
}

char* convert_iseq_to_cseq(int* iseq)
{
	int len = 0;
	while(iseq[len] > 0) {
		len++;
	}
	char* cseq = icalloc(len, char);
	for(int i = 0; i < len; i++)
	{
		int idx = iseq[i];
		if(idx >= NUM_AA_SCORE_MATRIX)
		{
			idx = NUM_AA_SCORE_MATRIX - 1;
		}
		cseq[i] = SCORE_MATRIX_AAS[idx];
	}
	cseq[len - 1 ] = '\0';
	return cseq;
}

int max3(int i, int j, int k)
{
	if(i>=j && i >=k)
	{

		return i;
	}
	else if (j >= k)
	{

		return j;
	}
	else{

		return k;
	}
}
// implement Needleman-Wunsch algorithm
// just get the score
int calculate_global_alignment_score(const char* seq1, const char* seq2)
{
	int* iseq1 = convert_cseq_to_iseq(seq1);
	int* iseq2 = convert_cseq_to_iseq(seq2);
	int len1 = strlen(seq1);
	int len2 = strlen(seq2);
	static int matrix[MAX_CDR3_LENGTH+1][MAX_CDR3_LENGTH+1];
	int gap = BLOSUM62_SCORE_MATRIX[0][23];

	for(int i = 0; i <= len1; i++)
	{
		matrix[0][i] = i*gap;
	}
	for(int i = 0; i <= len2; i++)
	{
		matrix[i][0] = i*gap;
	}
	int valueCross = 0;
	int valueLeft =0;
	int valueUp = 0;

	for(int i = 1; i <= len1; i++)
	{
		int idx_i = iseq1[i-1];
		for(int j = 1; j <= len2; j++)
		{
			int idx_j = iseq2[j-1];
			valueCross = matrix[j-1][i-1] + BLOSUM62_SCORE_MATRIX[idx_i][idx_j];
			valueLeft = matrix[j][i-1] + gap;
			valueUp = matrix[j-1][i] + gap;
			matrix[j][i] = max3(valueCross, valueLeft, valueUp);
		}
	}
	safe_free(iseq1);
	safe_free(iseq2);

	// debug("\n%s\n%s\n%d\n\n", seq1, seq2, matrix[len2][len1]);

	return matrix[len2][len1];
}

int get_distance(const char* s1, const char* s2)
{
	if(strcmp(s1, NA_SYMBOL) == 0)
	{
		if(strcmp(s2, NA_SYMBOL) == 0)
		{
			return NA_NA_DISTANCE;
		}else{
			return NA_CDR3_DISTANCE;
		}
	}else{
		if(strcmp(s2, NA_SYMBOL) == 0)
		{
			return NA_CDR3_DISTANCE;
		}else{
			return calculate_global_alignment_score(s1, s2);
		}
	}
}

void sort_cluster(si_vector *vec)
{
	int num_cdr3 = 0;
	for(int i = 0; i < kv_size(*vec); i++)
	{
		si_t *obj = kv_A(*vec, i);
		if (strcmp(obj->key, NA_SYMBOL) != 0)
		{
			num_cdr3++;
		}
	}
	if(num_cdr3 < 3)
	{
		sort_si_vector(vec);
		return;
	}

	i3_vector dist_vec;
	kv_init(dist_vec);

	for(int i = 0; i < kv_size(*vec); i++)
	{
		si_t *obj_i = (si_t*)kv_A(*vec, i);
		char* cdr3_i = obj_i->key;
		for(int j = i+1; j < kv_size(*vec); j++)
		{
			si_t *obj_j = (si_t*)kv_A(*vec, j);
			char* cdr3_j = obj_j->key;
			i3_t *obj = icalloc(1, i3_t);
			obj->idx1 = i;
			obj->idx2 = j;
			obj->dist = get_distance(cdr3_i, cdr3_j);
			kv_push(i3_t*, dist_vec, obj);
			// debug("i:%d j:%d\n", i, j);
		}
	}
	sort_i3_vector(&dist_vec);

	node_t* cur_node = 0;
	it_hash *dir_hash = kh_init(itm);
	i_vector ivec;
	// dbg_print("dist_vec size %d", kv_size(dist_vec));
	int tracking = 0;
	for(int i = 0; i < kv_size(dist_vec); i++)
	{
		i3_t* obj = (i3_t*)kv_A(dist_vec, i);
		// debug("idx1 %d idx2 %d dist %d\n", obj->idx1, obj->idx2, obj->dist);
		// debug("idx1:%d idx2:%d dist:%d\n", obj->idx1, obj->idx2, obj->dist);
		node_t* node1 = get_value_it_hash(dir_hash, obj->idx1);
		if(node1 == NULL)
		{
			node1 = create_node();
			tracking++;
			// create leaf node
			node1->idx = obj->idx1;
			// debug("====> node1->idx %d\n", obj->idx1);
		}
		node_t* node2 = get_value_it_hash(dir_hash, obj->idx2);

		if(node2 == NULL)
		{
			node2 = create_node();
			tracking++;
			// create leaf node
			node2->idx = obj->idx2;
			// debug(">>>>> node2->idx %d\n", obj->idx2);
		}
		// important logic hole
		if(node1 == node2) {
			continue;
		}

		cur_node = create_node();
		// create branch node
		cur_node->left = node1;
		cur_node->right = node2;

		kv_init(ivec);
		find_all_leaf_index(cur_node, &ivec);
		for(int j = 0; j < kv_size(ivec); j++)
		{
			int idx = kv_A(ivec, j);
			append_to_it_hash(dir_hash, idx, cur_node);
		}
		//clear ivec
		destruct_i_vector(&ivec);
	}
	destruct_it_hash(dir_hash);
	destruct_i3_vector(&dist_vec);
	kv_init(ivec);
	find_all_leaf_index(cur_node, &ivec);
	si_vector tmp; // holder all entries in the input vector
	kv_init(tmp);
	for(int i = 0; i < kv_size(*vec); i++)
	{
		si_t *obj = kv_A(*vec, i);
		kv_push(si_t*, tmp, obj);
		// debug("i: %d", obj->val);
	}
	// dbg_print("ivec size %d vec size %d tracking %d\n", kv_size(ivec), kv_size(*vec), tracking);

	for(int i = 0; i < kv_size(ivec); i++)
	{
		int idx = kv_A(ivec, i);
		si_t* obj = kv_A(tmp, idx);
		// debug("i: %d", obj->val);
		kv_A(*vec, i) = obj;
	}
	// clear ivec
	destruct_i_vector(&ivec);
	// clear tmp
	kv_destroy(tmp);
	// clear tree root and all children
	destruct_node(cur_node);
}
