#include <Python.h>
#include <stdbool.h>

const int EQ=2, MM=-4, GO=-6, GE=-2;

typedef struct DP_item{
    int score;
    bool from_ins;
    bool from_del;
} DPitem;

int str_len(char* s){
    int len=0;

    while (s[len]!='\0')
        len++;

    return len;
}

DPitem get_new_entry(DPitem left, DPitem left_up, DPitem up, bool match){

    int new_scores[3];

    if(left.from_del)
        new_scores[0] = left.score + GE;
    else
        new_scores[0] = left.score + GO;

    if(up.from_ins)
        new_scores[1] = up.score + GE;
    else
        new_scores[1] = up.score + GO;

    if(match)
        new_scores[2] = left_up.score + EQ;
    else
        new_scores[2] = left_up.score + MM;

    int i;
    DPitem new_item;
    new_item.score = 0;

    for(i=0;i<3;i++)
        if(new_scores[i]>new_item.score)
            new_item.score = new_scores[i];

    if(new_item.score==new_scores[0])
        new_item.from_del = true;
    else
        new_item.from_del = false;

    if(new_item.score==new_scores[1])
        new_item.from_ins = true;
    else
        new_item.from_ins = false;

    return new_item;

}

static PyObject *method_get_max_array(PyObject *self, PyObject *args){
    char *query, *ref;
    int i,j;

    if(!PyArg_ParseTuple(args, "ss", &query, &ref)){
        return NULL;
    }

    int query_len = str_len(query);
    int ref_len = str_len(ref);

    DPitem **matrix;
    matrix = (DPitem **)malloc(2*sizeof(DPitem *));

    for(i=0;i<2;i++)
        matrix[i] = (DPitem *)malloc((ref_len+1)*sizeof(DPitem));

    DPitem default_item;
    default_item.score = 0;
    default_item.from_del = false;
    default_item.from_ins = false;

    for(i=0;i<2;i++)
        for(j=0;j<ref_len+1;j++)
            matrix[i][j] = default_item;

    int current=1, previous=0;
    int *max_array = (int *)malloc(2*(query_len+1)*sizeof(int));
    max_array[0] = 0;
    max_array[1] = 0;

    for(i=0;i<query_len;i++){
        int max_score=0, max_index=0;
        for(j=0;j<ref_len;j++){
            bool match = query[i] == ref[j];
            matrix[current][j+1] = get_new_entry(matrix[current][j], matrix[previous][j], matrix[previous][j+1], match);

            if(matrix[current][j+1].score>max_score){
                max_score = matrix[current][j+1].score;
                max_index = j+1;
            }
        }

        max_array[2*(i+1)] = max_score;
        max_array[2*(i+1)+1] = max_index;

        current = (current+1)%2;
        previous = (previous+1)%2;
    }

    PyObject *result = PyList_New(query_len+1);
    for(i=0;i<query_len+1;i++){
        PyObject *item = PyList_New(2);
        PyObject *score = PyLong_FromLong(max_array[2*i]);
        PyObject *index = PyLong_FromLong(max_array[2*i+1]);

        PyList_SetItem(item, 0, score);
        PyList_SetItem(item, 1, index);

        PyList_SetItem(result, i, item);
    }

    // free
    for(i=0;i<2;i++) {
        free(matrix[i]);
    }
    free(matrix);
    free(max_array);

    return result;
}

static PyMethodDef get_max_array_Method[] = {
    {"get_max_array", method_get_max_array, METH_VARARGS, "get_max_array in C"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef testmodule = {
    PyModuleDef_HEAD_INIT,
    "dp_module",
    "C implementation for DP in SVanchor.",
    -1,
    get_max_array_Method
};

PyMODINIT_FUNC PyInit_dp_module(void){
    return PyModule_Create(&testmodule);
}
