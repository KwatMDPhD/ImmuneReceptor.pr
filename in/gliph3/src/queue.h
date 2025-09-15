#ifndef __QUEUE_H_
#define __QUEUE_H_
#include <stdbool.h>

#define MAXQUEUE 10000

typedef struct qnode
{
	int data;
	struct qnode* next;
}qnode_t;

typedef struct queue {
	qnode_t *front;
	qnode_t *rear;
	int items;
}queue_t;

void init_queue(queue_t *pq);
bool queue_is_full(const queue_t *pq);
bool queue_is_empty(const queue_t *pq);
int queue_item_count(const queue_t *pq);
bool enqueue(queue_t* pq, int data);
bool dequeue(queue_t* pq, int *data);
void empty_queue(queue_t* pq);
#endif
