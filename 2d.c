#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

//A point in 2d space
typedef struct{
    double x;
    double y;
} point;

//A node in our SBT
typedef struct __node__ {
    point pivot;
    int level;
    int leaves;
    struct __node__* parent;
    struct __node__* l;
    struct __node__* r;
    struct __node__* d;
} node;

int cmp_x(const void* x1, const void* x2){
    point p1 = *((point*)x1);
    point p2 = *((point*)x2);
    if (p1.x<p2.x) return -1;
    if (p1.x>p2.x) return  1;
    return 0;
}

int cmp_y(const void* x1, const void* x2){
    point p1 = *((point*)x1);
    point p2 = *((point*)x2);
    if (p1.y<p2.y) return -1;
    if (p1.y>p2.y) return  1;
    return 0;
}

int cmp(const void* x1, const void* x2){
    double y1 = *((double*)x1);
    double y2 = *((double*)x2);
    if (y1<y2) return -1;
    if (y1>y2) return  1;
    return 0;
}

void sort(double* xs, int len){
    qsort(xs,len,sizeof(double),cmp);
}

void sort_x(point* xs, int len){
    qsort(xs,len,sizeof(point),cmp_x);
}

void sort_y(point* xs, int len){
    qsort(xs,len,sizeof(point),cmp_y);
}

//Builds a y-dimension BST from the sub array ps_[start,end)
node* build_ytree(point* ps_, int start, int end){
    printf("Building ytree from %d to %d\n", start, end);
    int i;
    int len;
    int l = end-start;
    point* ps = malloc((l)*sizeof(point));

    memcpy(ps, &ps_[start], (l)*sizeof(point));

    //len is just l rounded up to the next power of 2
    if ((l & (l-1)) != 0) {
        len = 1<<((8*sizeof(int))-__builtin_clz(l));
    } else {
        len = l;
    }

    printf("len=%d, l=%d\n",len,l);

    node** nodes = malloc(len*sizeof(node*));
    int* indices = malloc(len*sizeof(int));
    node* new_node;

    sort_y(ps, l);

    for (i=0; i<l; i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->pivot = ps[i];
        nodes[i]->leaves = 1;
        indices[i] = i;
    }

    for (;i<len;i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->pivot.x = HUGE_VAL;
        nodes[i]->pivot.y = HUGE_VAL;
        nodes[i]->leaves = 0;
        indices[i] = i;
    }

    while (len>0){
        //printf("%d\n", len);
        len = len>>1;
        for (i=0; i<len; i++){
            new_node = malloc(sizeof(node));
            nodes[2*i]->parent   = new_node;
            nodes[2*i+1]->parent = new_node;

            new_node->l = nodes[2*i];
            new_node->r = nodes[2*i+1];
            new_node->d = NULL;
            new_node->leaves = new_node->l->leaves + new_node->r->leaves;
            new_node->pivot = ps[(indices[2*i]+indices[2*i+1])/2];
            //printf("len=%d, i=%d, indices[]...=%d\n", len, i, (indices[2*i]+indices[2*i+1])/2);
            indices[i] = (indices[2*i]+indices[2*i+1])/2;
            nodes[i] = new_node;
        }
    }

    free(indices);

    return nodes[0];
}

//Builds an x-dimension BST from the ps of length len
node* build_xtree(point* ps, int l){
    int i;
    int len;

    //len is just l rounded up to the next power of 2
    if ((l & (l-1)) != 0) {
        len = 1<<((8*sizeof(int))-__builtin_clz(l));
    } else {
        len = l;
    }

    node** nodes = malloc(len*sizeof(node*));
    int* indices = malloc(len*sizeof(int));
    int* starts = malloc(len*sizeof(int));
    int* ends = malloc(len*sizeof(int));
    node* new_node;


    sort_x(ps, l);

    for (i=0; i<l; i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->pivot = ps[i];
        nodes[i]->leaves = 1;
        indices[i] = i;
        starts[i] = i;
        ends[i] = i+1;
        nodes[i]->d = build_ytree(ps, i, i+1);
        printf("(%d,%d)\n", starts[i], ends[i]);
    }

    for (;i<len;i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->pivot.x = HUGE_VAL;
        nodes[i]->pivot.y = HUGE_VAL;
        nodes[i]->leaves = 0;
        indices[i] = i;
        starts[i] = l;
        ends[i] = l;
        printf("(%d,%d)\n", starts[i], ends[i]);
    }

    while (len>0){
        //printf("%d\n", len);
        len = len>>1;
        for (i=0; i<len; i++){
            new_node = malloc(sizeof(node));
            nodes[2*i]->parent   = new_node;
            nodes[2*i+1]->parent = new_node;

            new_node->l = nodes[2*i];
            new_node->r = nodes[2*i+1];
            new_node->d = build_ytree(ps, starts[2*i], ends[2*i+1]);
            new_node->leaves = new_node->l->leaves + new_node->r->leaves;
            new_node->pivot = ps[(indices[2*i]+indices[2*i+1])/2];
            //printf("len=%d, i=%d, indices[]...=%d\n", len, i, (indices[2*i]+indices[2*i+1])/2);
            printf("start=%d, end=%d\n", starts[2*i], ends[2*i+1]);
            starts[i] = starts[2*i];
            ends[i] = ends[2*i+1];
            indices[i] = (indices[2*i]+indices[2*i+1])/2;
            nodes[i] = new_node;
        }
    }

    free(indices);

    return nodes[0];
}

void print_ytree(node* n, int depth){
    if (n == NULL) return;
    for (int i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d\n", n->pivot.y, n->leaves);
    print_ytree(n->l, depth+1);
    print_ytree(n->r, depth+1);
}

void print_xtree(node* n, int depth){
    if (n == NULL) return;
    for (int i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d\n", n->pivot.x, n->leaves);
    for (int i = 0; i<depth; i++) printf(">");
    printf("ytree:\n");
    for (int i = 0; i<depth; i++) printf(">");
    printf("========\n");
    print_ytree(n->d, depth);
    for (int i = 0; i<depth; i++) printf(">");
    printf("========\n");
    print_xtree(n->l, depth+1);
    print_xtree(n->r, depth+1);
}


int main(void){
    double xs[15] = {4.0,14.0,2.0,11.0,0.0,1.0,7.0,10.0,15.0,3.0,6.0,5.0,12.0,8.0,9.0}; //Protip: x=13 is missing
    point* ps = malloc(15*sizeof(point));

    for (int i=0; i<15; i++){
        ps[i].x = xs[i];
        ps[i].y = (double) i;
    }

    node* root = build_xtree(ps, 15);

    print_xtree(root, 0);

    return 0;
}
