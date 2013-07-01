#include<stdio.h>
#include<stdlib.h>


//A node in our SBT
typedef struct __node__ {
    double pivot;
    int leaves;
    struct __node__* parent;
    struct __node__* l;
    struct __node__* r;
} node;

//A point in 2d space
typedef struct{
    double x;
    double y;
} point;

void test_recall(void){
    node* n = malloc(sizeof(node));
    int err = 0;

    printf("Testing recall...\n");

    n->pivot = 5;
    n->l = n;
    n->r = n;

    err &= (n->l->pivot == 5);
    err &= (n->r->pivot == 5);

    if (err == 0){
        printf("Yup, that's all fine\n");
    } else {
        printf("Nah, that's borked\n");
    }
}

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


//Builds a BST from the xs of length len assuming len is equal to a power of two
node* build_tree(double* xs, int len){
    node** nodes = malloc(len*sizeof(node*));
    int* indices = malloc(len*sizeof(int));
    node* new_node;

    sort(xs, len);

    for (int i=0; i<len; i++){
        nodes[i] = malloc(sizeof(node));
        nodes[i]->l = NULL;
        nodes[i]->r = NULL;
        nodes[i]->pivot = xs[i];
        nodes[i]->leaves = 1;
        indices[i] = i;
    }

    while (len>0){
        //printf("%d\n", len);
        len = len>>1;
        for (int i=0; i<len; i++){
            new_node = malloc(sizeof(node));
            nodes[2*i]->parent   = new_node;
            nodes[2*i+1]->parent = new_node;

            new_node->l = nodes[2*i];
            new_node->r = nodes[2*i+1];
            new_node->leaves = new_node->l->leaves + new_node->r->leaves;
            new_node->pivot = xs[(indices[2*i]+indices[2*i+1])/2];
            //printf("len=%d, i=%d, indices[]...=%d\n", len, i, (indices[2*i]+indices[2*i+1])/2);
            indices[i] = (indices[2*i]+indices[2*i+1])/2;
            nodes[i] = new_node;
        }
    }

    free(indices);

    return nodes[0];
}

void print_dft(node* n, int depth){
    if (n == NULL) return;
    for (int i = 0; i<depth; i++) printf(">");
    printf("%lf, leaves = %d\n", n->pivot,n->leaves);
    print_dft(n->l, depth+1);
    print_dft(n->r, depth+1);
}


int main(void){
    double xs[16] = {4.0,14.0,2.0,11.0,0.0,1.0,7.0,10.0,15.0,3.0,6.0,5.0,12.0,8.0,9.0,13.0};
    node* root = build_tree(xs, 16);

    print_dft(root, 0);

    return 0;
}
