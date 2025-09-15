#ifndef _HEAP_H_
#define _HEAP_H_

#include <stdio.h>
#include <stdlib.h>

typedef struct {
	float key;
	char* data;
	int index;
} hnode_t;

typedef struct {
	int size;
	hnode_t **elem;
} heap_t;


int parent(int i);

// to get index of left child of node at index i
int lchild(int i);

// to get index of right child of node at index i
int rchild(int i);

heap_t* init_heap();

void swap(hnode_t **n1, hnode_t **n2);

void heapify(heap_t* hp, int i);

void insert_heap_node(heap_t *hp, hnode_t* node);

hnode_t* pop_heap_node(heap_t* hp);

void destroy_heap(heap_t* hp);

void free_heap_node(hnode_t* node);

#endif
