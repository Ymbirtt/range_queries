#include<stdio.h>
#include<stdlib.h>

typedef struct __node__ {
    int pivot;
    struct __node__* l;
    struct __node__* r;
} node;

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

int main(void){
    double xs[10] = {5.0,3.0,1.0,2.0,8.0,4.0,7.0,6.0,9.0,10.0};

    sort(xs,10);

    for(int i = 0; i<10; i++){
        printf("%lf\n", xs[i]);
    }

    return 0;
}
