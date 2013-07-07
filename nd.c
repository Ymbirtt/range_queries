#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

typedef struct{
    double total_weight;
    int total_leaves;
} response;

typedef struct{
    double* x;
    int n;
} point;

//A node in our BST
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
    double s = 0.0;
    int i;
    for (i=0; i<p.n; i++){
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
    int i;
    int len;
    int l = end-start;
    point* ps = malloc((l)*sizeof(point));

    if (dim >= ps_[0].n) return NULL;
    printf("Building dimension %d tree from %d to %d\n", dim, start, end);


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
    int i;
    if (n == NULL) return;
    if (n->level >= max_level) return;


    printf("%d: ", n->level);
    for (i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d, weight = %.1lf\n", n->pivot, n->leaves, n->weight);
    if (n->d != NULL){
        printf("%d: ", n->level);
        for (i = 0; i<depth; i++) printf(">");
        printf("subtree:\n");
        printf("%d: ", n->level);
        for (i = 0; i<depth; i++) printf(">");
        printf("========\n");
        print_dtree(n->d, depth, max_level);
        printf("%d: ", n->level);
        for (i = 0; i<depth; i++) printf(">");
        printf("========\n");
    }
    print_dtree(n->l, depth+1, max_level);
    print_dtree(n->r, depth+1, max_level);
}

response query_layer(node* root, point lowers, point uppers, int layer){
    int i;
    node* u = root; //the node we reach when searching for a lower bound
    node* v = root; //the node we reach when searching for an upper
    response r;
    response r_;
    r.total_weight = 0;
    r.total_leaves = 0;

    for (i = 0; i<layer; i++) printf(">");
    printf("Querying layer %d\n", layer);

    //Navigate down the tree until upper and lower search paths diverge
    while (u == v){
        //If we are at a leaf node at this stage, interestingness happens
        if (u->l == NULL){
            //If this leaf's pivot is in the valid range, and we're at the lowest layer, then we can start adding things
            if (u->pivot > lowers.x[layer] && u->pivot < uppers.x[layer] && layer == lowers.n-1){
                for (i = 0; i<layer; i++) printf(">");
                printf("Hit leaf while finding z, adding value %.1lf\n", u->weight);
                r.total_leaves = u->leaves;
                r.total_weight = u->weight;
                return r;
            //If this leaf's pivot is in the valid range but we're not at the lowest layer, we can start recursing downwards
            } else if (u->pivot > lowers.x[layer] && u->pivot < uppers.x[layer]){
                //It's probably faster to do this iteratively, but it's much easier to do it recursively
                for (i = 0; i<layer; i++) printf(">");
                printf("Hit leaf while finding z\n");
            } else {
                return r;
            }
        }

        if (lowers.x[layer] < u->pivot){
            u = u->l;
            for (i = 0; i<layer; i++) printf(">");
            printf("u left, discriminator is %.1lf\n", u->pivot);
        } else {
            u = u->r;
            for (i = 0; i<layer; i++) printf(">");
            printf("u right, discriminator is %.1lf\n", u->pivot);
        }

        if (uppers.x[layer] < v->pivot){
            v = v->l;
            for (i = 0; i<layer; i++) printf(">");
            printf("v left\n");
        } else {
            v = v->r;
            for (i = 0; i<layer; i++) printf(">");
            printf("v right\n");
        }
    }

    for (i = 0; i<layer; i++) printf(">");
    printf("Found z\n");
    for (i = 0; i<layer; i++) printf(">");
    printf("Discriminators are %.1lf and %.1lf\n", u->pivot, v->pivot);
    for (i = 0; i<layer; i++) printf(">");
    printf("Underlying weights are %.1lf and %.1lf\n", u->weight, v->weight);

    if (layer == lowers.n-1){
        //If this is the bottom layer, add in the weights you find

        //Search for the lower bound - terminate at a leaf
        //Since construction always builds a perfectly balanced tree, we need only check one child
        //Whenever we go left, add in the total of the right child
        while (u->l != NULL){
            if (lowers.x[layer] < u->pivot){
                for (i = 0; i<layer; i++) printf(">");
                printf("u left\n");
                if (u->r != NULL){
                    r.total_weight += u->r->weight;
                    r.total_leaves += u->r->leaves;
                    for (i = 0; i<layer; i++) printf(">");
                    printf("Adding weight %.1lf\n", u->r->weight);
                }
                u = u->l;
            } else {
                u = u->r;
                for (i = 0; i<layer; i++) printf(">");
                printf("u right\n");
            }
        }
        //And don't forget that last point at the leaf
        if (u->pivot >= lowers.x[layer]){
            for (i = 0; i<layer; i++) printf(">");
            printf("Adding weight %.1lf\n", u->weight);
            r.total_weight += u->weight;
            r.total_leaves ++;
        }

        //Search for the upper bound - terminate at a leaf
        //Whenever we go right, add in the total of the left child
        while (v->l != NULL){
            if (uppers.x[layer] < v->pivot){
                v = v->l;
                for (i = 0; i<layer; i++) printf(">");
                printf("v left\n");
            } else {
                if (u->r != NULL){
                    r.total_weight += v->l->weight;
                    r.total_leaves += v->l->leaves;
                    for (i = 0; i<layer; i++) printf(">");
                    printf("Adding weight %.1lf\n", v->l->weight);
                }
                v = v->r;
                for (i = 0; i<layer; i++) printf(">");
                printf("v right\n");
            }
        }
        if (v->pivot < uppers.x[layer]){
            r.total_weight += v->weight;
            for (i = 0; i<layer; i++) printf(">");
            printf("Adding weight %.1lf\n", v->weight);
            r.total_leaves ++;
        }

        return r;
    } else {
        //If this is not the bottom layer, add in the weights of the subtrees below this one
        //
        //Do I need to store weights on layers other than the bottom? Probably not
        while (u->l != NULL){
            if (lowers.x[layer] < u->pivot){
                for (i = 0; i<layer; i++) printf(">");
                printf("u left\n");
                r_ = query_layer(u->r->d, lowers, uppers, layer+1);
                r.total_weight += r_.total_weight;
                r.total_leaves += r_.total_leaves;
                u = u->l;
            } else {
                u = u->r;
                for (i = 0; i<layer; i++) printf(">");
                printf("u right\n");
            }
        }
        //If this leaf node is valid in the relevant dimension, recurse and continue its consideration
        if(u->pivot > lowers.x[layer]){
            r_ = query_layer(u, lowers, uppers, layer+1);
            r.total_weight += r_.total_weight;
            r.total_leaves += r_.total_leaves;
        }
        for (i = 0; i<layer; i++) printf(">");
        printf("Found lower bound\n");
        while (v->l != NULL){
            if (uppers.x[layer] < v->pivot){
                v = v->l;
                for (i = 0; i<layer; i++) printf(">");
                printf("v left\n");
            } else {
                for (i = 0; i<layer; i++) printf(">");
                printf("v right\n");
                r_ = query_layer(v->l->d, lowers, uppers, layer+1);
                r.total_weight += r_.total_weight;
                r.total_leaves += r_.total_leaves;
                v = v->r;
            }
        }
        if(v->pivot < uppers.x[layer]){
            r_ = query_layer(v, lowers, uppers, layer+1);
            r.total_weight += r_.total_weight;
            r.total_leaves += r_.total_leaves;
        }
        for (i = 0; i<layer; i++) printf(">");
        printf("Found upper bound\n");

        return r;
    }
}

//Returns the average weight of all points strictly between lower and upper. Whilst I could and perhaps should make this
//inclusive, the fact that we'll be using this for continuous random numbers for integration means that this is low on
//the feature list
double query(node* root, point lower, point upper){
    response r = query_layer(root, lower, upper, 0);
    printf("Response found - total leaves = %d, total weight = %.1lf\n", r.total_leaves, r.total_weight);
    return r.total_weight/(double) r.total_leaves;
}

//Returns the same as above, but with a much simpler, much more obviously correct, but much slower method. Handy for
//testing
double stupid_query(point* ps, int len, point lower, point upper){
    int    total_points = 0;
    double total_weight = 0.0;
    int i,j;
    int valid = 1;

    //Loop through all the points
    for (i=0; i<len; i++){
        //Loop through all the dimensions until the number is invalid
        for (j=0; (j<lower.n) && valid; j++){
            if(!(ps[i].x[j] > lower.x[j] && ps[i].x[j] < upper.x[j])){
                valid = 0;
            }
        }
        if (valid){
            total_points++;
            total_weight += w(ps[i]);
        }
        valid = 1;
    }

    return total_weight/(double) total_points;
}

int main(void){
    double xs[15] = {4.0,14.0,2.0,11.0,0.0,1.0,7.0,10.0,15.0,3.0,6.0,5.0,12.0,8.0,9.0}; //Protip: x=13 is missing
    point* ps = malloc(15*sizeof(point));
    double avg;
    int i;

    for (i=0; i<15; i++){
        ps[i].x = malloc(2*sizeof(double));
        ps[i].x[0] = xs[i];
        ps[i].x[1] = (double) i;
        ps[i].n = 2;
    }

    node* root = build_dtree(ps, 0, 15, 0);

    print_dtree(root, 0, 2);

    printf("Computing average for points between ");
    print_point(ps[0]);
    printf(" and ");
    print_point(ps[14]);
    printf(".\n");

    avg = query(root, ps[0], ps[14]);
    printf("Clever average = %lf\n", avg);
    avg = stupid_query(ps, 15, ps[0], ps[14]);
    printf("Stupid average = %lf\n", avg);

    return 0;
}
