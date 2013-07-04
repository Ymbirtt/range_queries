#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

typedef struct{
    double* x;
    int n;
} point;

//A node in our SBT
typedef struct __node__ {
    double pivot;             //The disciminator value at this node; all values to the left are strictly smaller in the relevant dimension
    int level;                //The level of the tree - the dimension over which points are discriminated
    int leaves;               //The total number of leaves this node has as descendants
    double weight;            //The total weight of all leaves this node has as descendants
    struct __node__* parent;  //The parent of this node - currently defined but unused
    struct __node__* l;       //The left child
    struct __node__* r;       //The right child
    struct __node__* d;       //The tree on the level below this one, discriminating over a different dimension
} node;

//This needs to be a dumb value so that the program segfaults if you don't initialise it
//Also, for the love of god, learn C++ at some point...
int CMP_DIM = -1;

//Prints a point in n-d space but does not print a newline
void print_point(point p){
    int i;
    printf("(");
    for (i=0; i<p.n-1; i++){
        printf("%.1lf,", p.x[i]);
    }
    printf("%.1lf)", p.x[i]);
}

int cmp_dim(const void* x1, const void* x2){
    point p1 = *((point*)x1);
    point p2 = *((point*)x2);
    if (p1.x[CMP_DIM]<p2.x[CMP_DIM]) return -1;
    if (p1.x[CMP_DIM]>p2.x[CMP_DIM]) return  1;
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

//Sorts the xs of length len along dimension dim
void sort_dim(point* xs, int dim, int len){
    CMP_DIM = dim;
    qsort(xs,len,sizeof(point),cmp_dim);
}

//TODO: I need to look up function pointers again - eventually these will be passed into the constructor as arguments

//A weighting function
double w(point p){
    double s;
    for (int i=0; i<p.n; i++){
        s+=p.x[i];
    }
    return s;
}

//A semigroup operator
double add(double x1, double x2){
    return x1+x2;
}

//Builds a BST from the sub array ps_[start,end), discriminating by dimension dim
//Assumes that the ps_ are non-empty, and that they all exist in the same dimension
//PLEASE do not try and call this with points existing in multiple different spaces
node* build_dtree(point* ps_, int start, int end, int dim){
    if (dim >= ps_[0].n) return NULL;
    printf("Building dimension %d tree from %d to %d\n", dim, start, end);

    int i;
    int len;
    int l = end-start;
    point* ps = malloc((l)*sizeof(point));

    //Given the amount of in-place sorting happening here, making a local copy is probably advisable
    memcpy(ps, &ps_[start], (l)*sizeof(point));

    //len is just l rounded up to the next power of 2
    if ((l & (l-1)) != 0) {
        len = 1<<((8*sizeof(int))-__builtin_clz(l));
    } else {
        len = l;
    }

    node** nodes = malloc(len*sizeof(node*));
    int* indices = malloc(len*sizeof(int));
    node* new_node;
    int* starts = malloc(len*sizeof(int));
    int* ends = malloc(len*sizeof(int));

    sort_dim(ps, dim, l);

    for (i=0; i<l; i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->level = dim;
        nodes[i]->pivot = ps[i].x[dim];
        nodes[i]->leaves = 1;
        nodes[i]->weight = w(ps[i]);
        starts[i] = i;
        ends[i] = i+1;
        nodes[i]->d = build_dtree(ps, i, i+1, dim+1);
        indices[i] = i;
        //printf("(%d,%d)\n", starts[i], ends[i]);
    }

    for (;i<len;i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->level = dim;
        nodes[i]->pivot = HUGE_VAL;
        nodes[i]->leaves = 0;
        nodes[i]->weight = 0.0;
        starts[i] = l;
        ends[i] = l;
        nodes[i]->d = NULL;
        indices[i] = i;
        //printf("(%d,%d)\n", starts[i], ends[i]);
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
            printf("start=%d, end=%d\n", starts[2*i], ends[2*i+1]);
            new_node->d = build_dtree(ps, starts[2*i], ends[2*i+1],dim+1);
            new_node->level = dim;
            new_node->leaves = new_node->l->leaves + new_node->r->leaves;
            new_node->weight = add(new_node->l->weight, new_node->r->weight);
            new_node->pivot = ps[(indices[2*i]+indices[2*i+1])/2].x[dim];
            starts[i] = starts[2*i];
            ends[i] = ends[2*i+1];
            //printf("len=%d, i=%d, indices[]...=%d\n", len, i, (indices[2*i]+indices[2*i+1])/2);
            indices[i] = (indices[2*i]+indices[2*i+1])/2;
            nodes[i] = new_node;
        }
    }

    free(indices);

    return nodes[0];
}

void print_dtree(node* n, int depth, int max_level){
    if (n == NULL) return;
    if (n->level >= max_level) return;

    printf("%d: ", n->level);
    for (int i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d, weight = %.0lf\n", n->pivot, n->leaves, n->weight);
    if (n->d != NULL){
        printf("%d: ", n->level);
        for (int i = 0; i<depth; i++) printf(">");
        printf("subtree:\n");
        printf("%d: ", n->level);
        for (int i = 0; i<depth; i++) printf(">");
        printf("========\n");
        print_dtree(n->d, depth, max_level);
        printf("%d: ", n->level);
        for (int i = 0; i<depth; i++) printf(">");
        printf("========\n");
    }
    print_dtree(n->l, depth+1, max_level);
    print_dtree(n->r, depth+1, max_level);
}

void print_ytree(node* n, int depth){
    if (n == NULL) return;
    for (int i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d\n", n->pivot, n->leaves);
    print_ytree(n->l, depth+1);
    print_ytree(n->r, depth+1);
}

void print_xtree(node* n, int depth){
    if (n == NULL) return;
    //printf("%d:", n.level);
    for (int i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d, weight = %.0lf\n", n->pivot, n->leaves, n->weight);
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
        ps[i].x = malloc(2*sizeof(double));
        ps[i].x[0] = xs[i];
        ps[i].x[1] = (double) i;
        ps[i].n = 2;
    }

    node* root = build_dtree(ps, 0, 15, 0);

    print_dtree(root, 0, 2);

    return 0;
}
