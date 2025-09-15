#include "heap.h"
#include "ic_util.h"

int parent(int i) {
	return (i-1)/2;
}

// to get index of left child of node at index i
int lchild(int i) {
	return (2*i + 1);
}

// to get index of right child of node at index i
int rchild(int i) {
	return (2*i + 2);
}


heap_t* init_heap()
{
	heap_t* hp = malloc(sizeof(heap_t));
	hp->size = 0;
	hp->elem = 0;
	return hp;
}

void swap(hnode_t **n1, hnode_t **n2)
{
	hnode_t* temp = *n1;
	*n1 = *n2;
	*n2 = temp;
}

void heapify(heap_t* hp, int i)
{
	int smallest = i;
	if(lchild(i) < hp->size && hp->elem[lchild(i)]->key < hp->elem[i]->key)
	{
		smallest = lchild(i);
	}

	if(rchild(i) < hp->size && hp->elem[rchild(i)]->key < hp->elem[smallest]->key)
	{
		smallest = rchild(i);
	}

	if (smallest != i)
	{
		swap(&(hp->elem[i]), &(hp->elem[smallest]));
		heapify(hp, smallest);
	}
}

void insert_heap_node(heap_t *hp, hnode_t* node)
{
	if(hp->size)
	{
		hp->elem = realloc(hp->elem, (hp->size + 1) * sizeof(hnode_t*));
	}else{
		hp->elem = malloc(sizeof(hnode_t*));
	}
	int i = (hp->size)++;
	while (i && node->key < hp->elem[parent(i)]->key)
	{
		hp->elem[i] = hp->elem[parent(i)];
		i = parent(i);
	}
	hp->elem[i] = node;
}

hnode_t* pop_heap_node(heap_t* hp)
{
	if(hp->size)
	{
		hnode_t *node = hp->elem[0];
		hp->elem[0] = hp->elem[--(hp->size)];
		hp->elem = realloc(hp->elem, hp->size * sizeof(hnode_t*));
		heapify(hp, 0);
		return node;
	}
	else
	{
		return 0;
	}
}

void destroy_heap(heap_t* hp)
{
	for(int i = 0; i<hp->size; ++i)
	{
		free_heap_node(hp->elem[i]);
	}
	safe_free(hp->elem);
}

void free_heap_node(hnode_t* node)
{
	safe_free(node->data);
	safe_free(node);
}
