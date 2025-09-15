#include <stdio.h>
#include <stdlib.h>
#include "queue.h"


void init_queue(queue_t* pq)
{
	pq->front = pq->rear = NULL;
	pq->items = 0;
}

bool queue_is_full(const queue_t *pq)
{
	return pq->items == MAXQUEUE;
}

bool queue_is_empty(const queue_t* pq)
{
	return pq->items == 0;
}

int queue_item_count(const queue_t *pq)
{
	return pq->items;
}

bool enqueue(queue_t *pq, int data)
{
	qnode_t* pnew;
	if(queue_is_full(pq))
	{
		return false;
	}
	pnew = (qnode_t*) malloc(sizeof(qnode_t));
	if(pnew == NULL)
	{
		fprintf(stderr, "Unable to allocate memory!\n");
		exit(1);
	}
	pnew->data = data;
	pnew->next = NULL;
	if(queue_is_empty(pq))
	{
		pq->front = pnew;
	}else{
		pq->rear->next = pnew;
	}
	pq->rear = pnew;
	pq->items++;
	return true;
}

bool dequeue(queue_t * pq, int *data)
{
	qnode_t* pt;
	if(queue_is_empty(pq))
	{
		return false;
	}
	*data = pq->front->data;
	pt = pq->front;
	pq->front = pq->front->next;
	free(pt);
	pq->items--;
	if(pq->items == 0)
	{
		pq->rear = NULL;
	}
	return true;
}

void empty_queue(queue_t *pq)
{
	int dummy;
	while(!queue_is_empty(pq))
	{
		dequeue(pq, &dummy);
	}
}
