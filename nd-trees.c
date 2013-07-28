#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<time.h>

#include"mt.h"

typedef struct{
    double total_weight;
    int total_leaves;
} response;

typedef struct __point__{
    double* x;
    int n;
} point;


//Defines a box as a point on a corner, and d orthogonal vectors
//representing each of its edges
typedef struct{
    point corner;
    point* side;
} box;

//A node in our BST
typedef struct __node__ {
    double pivot;             //The disciminator value at this node; all values to the left are strictly smaller in the relevant dimension
    int dim;                  //The dimension over which points are discriminated
    int leaves;               //The total number of leaves this node has as descendants
    double weight;            //The total weight of all leaves this node has as descendants
    point p;                  //If this is a leaf node, then p will be the point at this leaf
    struct __node__* parent;  //The parent of this node - currently defined but unused
    struct __node__* l;       //The left child
    struct __node__* r;       //The right child
} node;

//This needs to be a dumb value so that the program segfaults if you don't initialise it
//Also, for the love of god, learn C++ at some point...
int CMP_DIM = -1;

//Prints a point in n-d space but does not print a newline
void print_point(point p){
    int i;
    if (p.n == -1) return;
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

//Sorts the ps of length len along dimension dim
void sort_dim(point* ps, int dim, int len){
    CMP_DIM = dim;
    qsort(ps,len,sizeof(point),cmp_dim);
}

//Sorts ps[start,end) in dimension dim
void sort_dim_range(point* ps, int start, int end, int dim){
    CMP_DIM = dim;
    qsort(&ps[start], end-start, sizeof(point), cmp_dim);
}

//Returns a continuous, uniormly distributed random-ish number between a and b
//seeding function should be called beforehand
//mmmm, delicious mersenne twister
double uniform(double a, double b){
    return genrand_real1()*(b-a) + a;
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

//Builds a BST from the sub array ps_[start,end), discriminating by dimension dim
//Assumes that the ps are non-empty, and that they all exist in the same dimension
//PLEASE do not try and call this with points existing in multiple different spaces
node* build_dtree(point* ps, int start, int end, int dim){
    node* new_node = malloc(sizeof(node));

    //printf("Building dimension %d tree from %d to %d\n", dim, start, end);

    if (start == end-1){
        //Leaf node
        new_node->leaves = 1;
        new_node->weight = w(ps[start]);
        new_node->l = NULL;
        new_node->r = NULL;
        new_node->p = ps[start];
    } else if (start == end) {
        //Empty node
        new_node->leaves = 0;
        new_node->weight = 0;
        new_node->l = NULL;
        new_node->r = NULL;
        new_node->p.n = -1;
    } else {
        //Inner node
        sort_dim_range(ps, start, end, dim);

        new_node->pivot = ps[(start+end)/2].x[dim];
        new_node->dim = dim;
        new_node->p.x = NULL;

        new_node->l = build_dtree(ps, start, (start+end)/2, (dim+1) % ps[0].n);
        new_node->r = build_dtree(ps, (start+end)/2, end, (dim+1) % ps[0].n);

        new_node->leaves = new_node->l->leaves + new_node->r->leaves;
        new_node->weight = new_node->l->weight + new_node->r->weight;

        new_node->l->parent = new_node;
        new_node->r->parent = new_node;
    }

    return new_node;
}

void print_dtree(node* n, int depth){
    int i;
    if (n == NULL) return;

    for (i = 0; i<depth; i++) printf(">");
    if (n->p.x == NULL){
        printf("%.1lf, split = %d, leaves = %d, weight = %.1lf\n", n->pivot, n->dim, n->leaves, n->weight);
    } else {
        print_point(n->p);
        printf(", weight = %.1lf\n", n->weight);
    }
    print_dtree(n->l, depth+1);
    print_dtree(n->r, depth+1);
}

//The cell coordinates specify the corners of the box that is currently under consideration, and is assumed to not be a
//subset of the query box
response query_cell(node* root, point lower, point upper, point cell_lower, point cell_upper){
    response r;
    response r_;
    int valid = 1;
    int i;
    double old_val;

    r.total_leaves = 0;
    r.total_weight = 0.0;

    if (root->leaves == 0){
        //Empty node
        return r;
    } else if (root->leaves ==1){
        //Leaf node

        //Is this point valid?
        for (i = 0; i<root->p.n && valid; i++){
            if (!(root->p.x[i] < upper.x[i] && root->p.x[i] > lower.x[i])){
                valid = 0;
            }
        }
        //If so, return it and its weight
        if (valid){
            r.total_leaves = 1;
            r.total_weight = root->weight;
        }
        return r;
    } else {
        //Internal node

        //Is the search cell a subset of the query cell?
        for (i = 0; i<lower.n && valid; i++){
            if (!(cell_upper.x[i] < upper.x[i] && cell_lower.x[i] > lower.x[i])){
                valid = 0;
            }
        }
        //If so, return the total weight of all children
        if (valid){
            r.total_leaves = root->leaves;
            r.total_weight = root->weight;
            return r;
        }
        //If not, start searching

        //Search the left subtree, this search cell is upper bounded in the relevant dimension by the pivot
        old_val = cell_upper.x[root->dim];
        cell_upper.x[root->dim] = root->pivot;
        r_ = query_cell(root->l, lower, upper, cell_lower, cell_upper);
        cell_upper.x[root->dim] = old_val;

        r.total_leaves += r_.total_leaves;
        r.total_weight += r_.total_weight;

        //Search the right subtree, this search cell is lower bounded in the relevant dimension by the pivot
        old_val = cell_lower.x[root->dim];
        cell_lower.x[root->dim] = root->pivot;
        r_ = query_cell(root->r, lower, upper, cell_lower, cell_upper);
        cell_lower.x[root->dim] = old_val;

        r.total_leaves += r_.total_leaves;
        r.total_weight += r_.total_weight;

        return r;
    }
}


//Returns the average weight of all points above lower and strictly below upper. Whilst I could and perhaps should make this
//inclusive, the fact that we'll be using this for continuous random numbers for integration means that this is low on
//the feature list
double query(node* root, point lower, point upper){
    response r;
    int i;
    point cell_lower, cell_upper;

    assert(lower.n == upper.n);

    cell_lower.n = lower.n;
    cell_upper.n = lower.n;

    cell_lower.x = malloc(lower.n*sizeof(double));
    cell_upper.x = malloc(lower.n*sizeof(double));

    for (i=0; i<lower.n; i++){
        //The search range is initially the entire space. We'll be tightening it.
        cell_lower.x[i] = -HUGE_VAL;
        cell_upper.x[i] =  HUGE_VAL;
    }

    r = query_cell(root, lower, upper, cell_lower, cell_upper);

    printf("Total points found: %d, Total weights: %lf\n", r.total_leaves, r.total_weight);
    return r.total_weight/(double) r.total_leaves;
}


//Returns the dot product between two vectors
//If this is positive, the angle between them is acute
//If negative, the angle is obtuse
double dot_prod(point v1, point v2){
    int i;
    double prod = 0.0;
    for (i=0; i<v1.n; i++){
        prod += v1.x[i] * v2.x[i];
    }
    return prod;
}

//Stores p1-p2 in dest
void point_subtract(point dest, point p1, point p2){
    int i;
    for (i=0; i<p1.n; i++){
        dest.x[i] = p1.x[i] - p2.x[i];
    }
}

//Stores p1+p2 in dest
void point_add(point dest, point p1, point p2){
    int i;
    for (i=0; i<p1.n; i++){
        dest.x[i] = p1.x[i] + p2.x[i];
    }
}

//Returns 0 iff the given point does not lie strictly within the given box
int point_inside(box b, point p){
    point t1, t2;
    int i;

    t1.x = malloc(p.n * sizeof(double));
    t1.n = p.n;

    t2.x = malloc(p.n * sizeof(double));
    t2.n = p.n;


    point_subtract(t1, p, b.corner);
    for (i=0; i<p.n; i++){
        if(dot_prod(t1, b.side[i]) <= 0.0){
            free(t1.x);
            free(t2.x);
            return 0;
        }
    }

    for (i=0; i<p.n; i++){
        point_add(t2, b.corner, b.side[i]);
        point_subtract(t1, p, t2);
        if(dot_prod(t1, b.side[i]) > 0.0){
            free(t1.x);
            free(t2.x);
            return 0;
        }
    }
    free(t1.x);
    free(t2.x);
    return 1;
}

//Worst.
//Method.
//Ever.
//Check if all 2^d points that define the given cell are inside b.
//If all of them are, then return -1
//If none of them are, return 1
//Otherwise, return 0
int cell_cmp(const box b, const point cell_lower, const point cell_upper){
    int i,corner,num_corners;
    point tmp;
    int all_outside = 1;
    int all_inside = 1;

    tmp.n = cell_lower.n;
    tmp.x = malloc(tmp.n * sizeof(double));
    
    num_corners = (1<<cell_lower.n);
    
    //Encode each corner of the cell as a bit-string
    for(corner = 0; corner < num_corners; corner++){
        //printf("%d\n", corner);
        for(i=0; i<cell_lower.n; i++){
            if((corner & (1<<i)) != 0){
                //bit i of corner is set
                tmp.x[i] = cell_upper.x[i];
            }else{
                //bit i of corner is unset
                tmp.x[i] = cell_lower.x[i];
            }
        }
        if (point_inside(b,tmp)){
            //print_point(tmp);
            all_outside = 0;
        } else {
            //print_point(tmp);
            all_inside = 0;
        }
        if (all_inside == 0 && all_outside == 0){
            free(tmp.x);
            return 0;
        }
    }
    
    if(all_inside){
        free(tmp.x);    
        return -1;
    } else {
        free(tmp.x);
        return 1;
    }
}

response query_box_cell(node* root, box b, point cell_lower, point cell_upper){
    response r;
    response r_;
    int valid = 1;
    double old_val;

    r.total_leaves = 0;
    r.total_weight = 0.0;

    if (root->leaves == 0){
        //Empty node
        return r;
    } else if (root->leaves ==1){
        //Leaf node
        //Is this point valid?
        //If so, return it and its weight
        if (point_inside(b, root->p)){
            r.total_leaves = 1;
            r.total_weight = root->weight;
        }
        return r;
    } else {
        //Internal node

        //compare the current cell to the box
        valid = cell_cmp(b, cell_lower, cell_upper);

        //If the cell is completely inside the box, return the total weight of
        //all its children
        if (valid < 0){
            r.total_leaves = root->leaves;
            r.total_weight = root->weight;
            return r;
        }

        //If the cell is completely outside of the box, return 0
        if (valid > 0){
            return r;
        }

        //Otherwise, continue searching

        //Search the left subtree, this search cell is upper bounded in the relevant dimension by the pivot
        old_val = cell_upper.x[root->dim];
        cell_upper.x[root->dim] = root->pivot;
        r_ = query_box_cell(root->l, b, cell_lower, cell_upper);
        cell_upper.x[root->dim] = old_val;

        r.total_leaves += r_.total_leaves;
        r.total_weight += r_.total_weight;

        //Search the right subtree, this search cell is lower bounded in the relevant dimension by the pivot
        old_val = cell_lower.x[root->dim];
        cell_lower.x[root->dim] = root->pivot;
        r_ = query_box_cell(root->r, b, cell_lower, cell_upper);
        cell_lower.x[root->dim] = old_val;

        r.total_leaves += r_.total_leaves;
        r.total_weight += r_.total_weight;

        return r;
    }
}



//Generates a random d-dimensional vector
point rand_vector(int d, double min, double max){
    int i;
    point p;

    p.n = d;
    p.x = malloc(d*sizeof(double));

    for (i=0; i<d; i++){
        p.x[i] = uniform(min,max);
    }

    return p;
}


//Stores the projection operator applied to (u,v) in dest
void proj(point dest, point u, point v){
    int i;
    double t1, t2;

    t1 = dot_prod(u,v);
    t2 = dot_prod(u,u);

    t1 = t1/t2;

    for(i=0;i<u.n;i++){
        dest.x[i] = u.x[i] * t1;
    }
}

//Generates a random d-dimensional box, with vectors uniformly distributed in
//(min, max)
box rand_box(int d, double min, double max){
    int i,j;
    box b;
    point* temp_points;
    point temp_point;

    b.corner = rand_vector(d, min, max);
    b.side = malloc(d*sizeof(point));
    temp_points = malloc(d*sizeof(point));
    temp_point.n = d;
    temp_point.x = malloc(d*sizeof(point));

    for (i=0; i<d; i++){
        b.side[i] = rand_vector(d, min, max);
        temp_points[i].n = d;
        temp_points[i].x = malloc(d*sizeof(double));
        for (j=0; j<d; j++){
            temp_points[i].x[j] = 0.0;
        }
    }

    for(i=0; i<d; i++){
        point_add(temp_points[i], b.side[i], temp_points[i]);
        for(j=0; j<i-1; j++){
            proj(temp_point, temp_points[j], b.side[i]);
            point_subtract(temp_points[i], temp_points[i], temp_point);
        }
    }

    for(i=0; i<d; i++){
        free(b.side[i].x);
    }
    
    free(b.side);
    free(temp_point.x);
    b.side = temp_points;

    return b;
}

//Returns the average weight of all points inside the box defined by b.
double query_box(node* root, box b){
    response r;
    int i;
    point cell_lower, cell_upper;

    cell_lower.n = b.corner.n;
    cell_upper.n = b.corner.n;

    cell_lower.x = malloc(b.corner.n*sizeof(double));
    cell_upper.x = malloc(b.corner.n*sizeof(double));

    for (i=0; i<b.corner.n; i++){
        cell_lower.x[i] = -HUGE_VAL;
        cell_upper.x[i] =  HUGE_VAL;
    }

    r = query_box_cell(root, b, cell_lower, cell_upper);
    
    free(cell_lower.x);
    free(cell_upper.x);
    
    //printf("Total points found: %d, Total weights: %lf\n", r.total_leaves, r.total_weight);

    return r.total_weight/(double) r.total_leaves;

}

double stupid_query_box(point* ps, int len, box b){
    int total_points = 0;
    double total_weight = 0.0;
    int i;

    for (i=0; i<len; i++){
        if (point_inside(b,ps[i])){
            total_points++;
            total_weight += w(ps[i]);
        }
    }
    //printf("Total points found: %d, Total weights: %lf\n", total_points, total_weight);
    return total_weight/(double) total_points;
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
            //print_point(ps[i]);
            //printf(" is valid \n");
        }
        valid = 1;
    }
    printf("Total points found: %d, Total weights: %lf\n", total_points, total_weight);
    return total_weight/(double) total_points;
}
int main(int argc, char** argv){
    
    if (argc != 2){
        printf("NOPE NOPE NOPE");
        return 0;
    }
    
    int dims = 2;
    int points = 1<<atoi(argv[1]);
    int queries = 1000;
    
    printf("%d,", points);
    
    point* ps = malloc(points*sizeof(point));

    point lowers;
    point uppers;

    box b;

    clock_t tic1 = 0;
    clock_t toc1 = 0;
    clock_t tic2 = 0;
    clock_t toc2 = 0;

    //double avg1, avg2;
    //double r1, r2;
    int i,j;

    node* root;

    init_genrand(0);
    //printf("Generating points\n");

    for (i=0; i<points; i++){
        ps[i].x = malloc(dims*sizeof(double));
        ps[i].n = dims;
        for (j=0; j<dims; j++){
            ps[i].x[j] = uniform(0,10);
        }
    }
    tic1 = clock();
    root = build_dtree(ps, 0, points, 0);
    toc1 += clock() - tic1;

    printf("%d,", toc1);

    toc1 = 0;

    //sort_dim(ps, 0, points);
    //printf("Points are:\n");
    //for (i=0; i<points; i++){
    //    print_point(ps[i]);
    //    printf("\n");
    //}

    //print_dtree(root, 0);

    //printf("Computing average for points between ");
    //print_point(lowers);
    //printf(" and ");
    //print_point(uppers);
    //printf(".\n");

    lowers.x = malloc(dims*sizeof(double));
    uppers.x = malloc(dims*sizeof(double));

    lowers.n = dims;
    uppers.n = dims;

    //printf("Performing %d searches\n", queries);

    for (i=0; i<queries; i++){
        //if(i%(queries/100) == 0) printf("%d\n",i);
        b = rand_box(dims, 0, 5);
        tic1 = clock();
        query_box(root, b);
        toc1 += clock() - tic1;

        tic2 = clock();
        stupid_query_box(ps, points, b);
        toc2 += clock()-tic2;
    }
    printf("%d,", toc1);
    printf("%d\n",toc2);

    return 0;
}
