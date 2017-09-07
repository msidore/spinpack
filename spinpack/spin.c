#include <stdio.h>
#include <stdlib.h>
#include <python2.7/Python.h>
#include <math.h>
#define PY_ARRAY_UNIQUE_SYMBOL get_spinangle_traj_ARRAY_API
/* This one makes the compilation crash if not commented, despite the warning message if it is */
//~ #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
/* Not typecasting anything */
#include "spin.h"

/* Docstring */
static char module_docstring[] =
    "This module provides an interface for calculating a spinangle using C.";
static char get_spinangle_traj_docstring[] =
    "Calculate the spinangle of some coordinates relative to reference coordinates.";

static PyObject *get_spinangle_traj(PyObject *self, PyObject *args){

    /* Doing it again - will take as input 2 matrix, one as the list of coordinates of one frame of the trajectory
     * and the second as the list of coordinates from the reference frame */

    /* So, the first one will be the matrix of trajectory coordinates */
    PyObject *coords_matrix;
    /* And the second one is the matrix of the reference coordinates */
    PyObject *ref_matrix;
    /* How many atoms ? */
    int num_atoms;
    /* debug */
    double debug=0.0;

    /* Should get the python objects */
    if (!PyArg_ParseTuple(args,"O!O!i",&PyArray_Type,&coords_matrix,&PyArray_Type,&ref_matrix,&num_atoms)){
        PyErr_SetString(PyExc_ValueError,"Error while parsing the trajectory coordinates in get_spinangle_traj");
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    double *coords    = (double*)PyArray_DATA(coords_matrix);
    double *ref    = (double*)PyArray_DATA(ref_matrix);

    /* Call the external C function to get the spinangle */
    double spinangle = get_spinangle(coords, ref, num_atoms, debug);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("f", spinangle);

    return ret;
}

double get_spinangle(double *coords, double *ref, int num_atoms, double debug){
    /* Do stuff now, and hope everything went smoothly */

    /* Initialize the (future) angle and the quaternion matrix */
    double angle = 0.0;
    int i, j;
    /* And malloc */
    double **quaternion_matrix = malloc(4*sizeof(double*));
    for (i=0; i<4; i++) {
        quaternion_matrix[i] = malloc(4*sizeof(double));
    }

    /* Build the quaternion matrix */
    build_matrix(coords, ref, quaternion_matrix, num_atoms, debug);

    /* Eigenvector and Eigenvalue matrix */
    double *eigval_matrix = malloc(4*sizeof(double));
    double **eigvec_matrix = malloc(4*sizeof(double*));
    for (i=0; i<4; i++) {
        eigvec_matrix[i] = malloc(4*sizeof(double));
    }
    /* Fill them with 0s */
    for (i=0; i<4; i++) {
        eigval_matrix[i] = 0.0;
        for (j=0; j<4; j++){
            eigvec_matrix[i][j] = 0.0;
        }
    }

    /* Diagnoalize the quaternion matrix */
    int jacobi_numrotations = 0;
    jacobi(quaternion_matrix, eigval_matrix, eigvec_matrix, num_atoms, jacobi_numrotations);

    /* Free the quaternion matrix */
    for (i=0; i<4; i++) {
        free(quaternion_matrix[i]);
    }
    free(quaternion_matrix);

    /* Sort the result */
    eigsrt(eigval_matrix, eigvec_matrix);

    /* Transpose the eigenvectors */
    transpose(eigvec_matrix);

    /* Normalize the eigenvectors */
    normalize(eigvec_matrix);

    /* Now we have the maximum eigenvalue and the maximum eigenvector */
    
    /* Maximum eigenvalue
       aka lambda
       commented because we actually don't use lambda */
    //~ double max_eigval = 0.0;
    //~ max_eigval = eigval_matrix[0];
    
    /* Maximum eigenvector
       aka q, the rotation itself, from the first row of eigenvector */
    double *max_eigvec = malloc(4*sizeof(double));
    for (i=0; i<4; i++){
        max_eigvec[i] = eigvec_matrix[0][i];
    }

    /* Now we only have to get the angle from this, around the default  (0, 0, 1) axis */
    get_angle_quaternion(max_eigvec, &angle);

    /* Free eigvec and eigval */
    free(eigval_matrix);
    for (i=0; i<4; i++) {
        free(eigvec_matrix[i]);
    }
    free(eigvec_matrix);
    free(max_eigvec);

    return angle;
}

void build_matrix(double *coords, double *ref, double** quaternion_matrix, int num_atoms, double debug){

    /* Build the quaternion matrix - first with the correlation matrix, with memory allocation */
    int i, j;
    double **cor_matrix = malloc(3*sizeof(double*));
    for (i=0; i<3; i++){
        cor_matrix[i] = malloc(3*sizeof(double));
    }

    /* Initialize the correlation matrix */
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            cor_matrix[i][j] = 0.0;
        }
    }

    /* Fill the correlation matrix */
    for (i=0; i<num_atoms*3; i+=3){
        cor_matrix[0][0] += coords[i] * ref[i];
        cor_matrix[0][1] += coords[i] * ref[i+1];
        cor_matrix[0][2] += coords[i] * ref[i+2];
        cor_matrix[1][0] += coords[i+1] * ref[i];
        cor_matrix[1][1] += coords[i+1] * ref[i+1];
        cor_matrix[1][2] += coords[i+1] * ref[i+2];
        cor_matrix[2][0] += coords[i+2] * ref[i];
        cor_matrix[2][1] += coords[i+2] * ref[i+1];
        cor_matrix[2][2] += coords[i+2] * ref[i+2];
    }

    /* And now the quaternion matrix */
    /* Initialize it */
    for (i=0; i<4; i++){
        for (j=0; j<4; j++){
            quaternion_matrix[i][j] = 0.0;
        }
    }

    /* And fill it */
    quaternion_matrix[0][0] = cor_matrix[0][0] + cor_matrix[1][1] + cor_matrix[2][2];
    quaternion_matrix[1][0] = cor_matrix[1][2] - cor_matrix[2][1];
    quaternion_matrix[0][1] = quaternion_matrix[1][0];
    quaternion_matrix[2][0] = cor_matrix[2][0] - cor_matrix[0][2];
    quaternion_matrix[0][2] = quaternion_matrix[2][0];
    quaternion_matrix[3][0] = cor_matrix[0][1] - cor_matrix[1][0];
    quaternion_matrix[0][3] = quaternion_matrix[3][0];
    quaternion_matrix[1][1] = cor_matrix[0][0] - cor_matrix[1][1] - cor_matrix[2][2];
    quaternion_matrix[2][1] = cor_matrix[0][1] + cor_matrix[1][0];
    quaternion_matrix[1][2] = quaternion_matrix[2][1];
    quaternion_matrix[2][2] = - cor_matrix[0][0] + cor_matrix[1][1] - cor_matrix[2][2];
    quaternion_matrix[3][1] = cor_matrix[0][2] + cor_matrix[2][0];
    quaternion_matrix[1][3] = quaternion_matrix[3][1];
    quaternion_matrix[3][2] = cor_matrix[1][2] + cor_matrix[2][1];
    quaternion_matrix[2][3] = quaternion_matrix[3][2];
    quaternion_matrix[3][3] = - cor_matrix[0][0] - cor_matrix[1][1] + cor_matrix[2][2];

    /* Return the quaternion matrix */
    return;
}

void rotate_stuff(double** a, int i, int j, int k, int l, double* g, double* h, double tau, double s){

    /* Function to rotate things in the matrix */

    *g = a[i][j];
    *h = a[k][l];
    a[i][j] = *g - s*(*h + *g*tau);
    a[k][l] = *h + s*(*g - *h*tau);

    return;

}

void jacobi(double** quaternion_matrix, double* eigval_matrix, double** eigvec_matrix, int num_atoms, int numb_rot){

    /* Jacobi for matrix diagonalization, row-wise */

    /* Start by initializing stuff */
    double *b = malloc(4*sizeof(double));
    double *z = malloc(4*sizeof(double));
    int i = 0.0, j = 0.0, i2 = 0.0, j2 = 0.0;
    /* Fill with zeroes */
    for (i=0; i<4; i++){
        b[i] = 0.0;
        z[i] = 0.0;
    }
    double sum = 0.0, treshold = 0.0, g = 0.0, h = 0.0, t = 0.0, theta = 0.0, c = 0.0, s = 0.0, tau = 0.0;

    /* Fill eigval and eigvec */
    for (i=0; i<4; i++){
        for (j=0; j<4; j++){
            eigvec_matrix[i][j] = 0.0;
        eigvec_matrix[i][i] = 1.0;
        }
    }
    /* Do it right by filling z with zeroes and b with what has to */
    for (i=0; i<4; i++){
        b[i] = eigval_matrix[i] = quaternion_matrix[i][i];
        z[i] = 0.0;
    }
    /* Start the big loop ! For each iteration, annihilate a member of the quaternion matrix */
    for (i2=0; i2<50; i2++){
        /* Get the sum of everything that's not in the diagonal */
        sum = 0.0;
        for (i=0; i<3; i++){
            for (j=i+1; j<4; j++){
                sum += fabs(quaternion_matrix[i][j]);
            }
        }
        /* If the sum doesn't increase, we're set, return ! */
        if (sum == 0.0){
            return;
        }
        /* Set the treshold */
        if (i2 < 4){
            treshold = 0.2*sum / (4*4);
        } else {
            treshold = 0.0;
        }
        /* Here the second big loop */
        for (i=0; i<3; i++){
            for (j=i+1; j<4; j++){
                g = 100.0 * fabs(quaternion_matrix[i][j]);
                if ( ( i2 > 4 ) && ( (fabs(eigval_matrix[i]) + g) == fabs(eigval_matrix[i]) ) && ( (fabs(eigval_matrix[j]) + g) == fabs(eigval_matrix[j]) )){
                    quaternion_matrix[i][j] = 0.0;
                } else if (fabs(quaternion_matrix[i][j]) > treshold) {
                    h = eigval_matrix[j] - eigval_matrix[i];
                    if (fabs(h) + g == fabs(h)) {
                        t = (quaternion_matrix[i][j]) / h;
                    } else {
                        theta = 0.5 * h / (quaternion_matrix[i][j]);
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
                        if (theta < 0.0){
                            t = -t;
                        }
                    }
                    c = 1.0 / sqrt(1 + t*t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t*quaternion_matrix[i][j];
                    z[i] -= h;
                    z[j] += h;
                    eigval_matrix[i] -= h;
                    eigval_matrix[j] += h;
                    /* There's the annihilation */
                    quaternion_matrix[i][j] = 0.0;
                    /* And the rotations, before ending the iteration */
                    for (j2=0; j2<i-1; j2++){
                        rotate_stuff(quaternion_matrix, j2, i, j2, j, &g, &h, tau, s);
                    }
                    for (j2=i+1; j2<j; j2++){
                        rotate_stuff(quaternion_matrix, i, j2, j2, j, &g, &h, tau, s);
                    }
                    for (j2=j+1; j2<4; j2++){
                        rotate_stuff(quaternion_matrix, i, j2, j, j2, &g, &h, tau, s);
                    }
                    for (j2=0; j2<4; j2++){
                        rotate_stuff(eigvec_matrix, j2, i, j2, j, &g, &h, tau, s);
                    }
                    numb_rot += 1;
                }
            }
        }
        for (i=0; i<4; i++){
            b[i] += z[i];
            eigval_matrix[i] = b[i];
            z[i] = 0.0;
        }
    }

    free(b);
    free(z);

    return;

}

void eigsrt(double* eigval_matrix, double** eigvec_matrix){

    /* Sort the eigenvectors and values */

    int i, j, k = 0;
    double p = 0.0;

    for (i=0; i<4; i++){
        p = eigval_matrix[i];
        k = i;
        for (j=i+1; j<4; j++){
            if (eigval_matrix[j] >= p){
                p = eigval_matrix[j];
                k = j;
            }
        }
        if (k != i){
            eigval_matrix[k] = eigval_matrix[i];
            eigval_matrix[i] = p;
            for (j=0; j<4; j++){
                p = eigvec_matrix[j][i];
                eigvec_matrix[j][i] = eigvec_matrix[j][k];
                eigvec_matrix[j][k] = p;
            }
        }
    }

    return;

}

void transpose(double** eigvec_matrix){

    /* Small routine to transpose a 4 4 matrix */

    int i, j;
    double p = 0.0;

    for (i=0; i<4; i++){
        for (j=i+1; j<4; j++){
            p = eigvec_matrix[i][j];
            eigvec_matrix[i][j] = eigvec_matrix[j][i];
            eigvec_matrix[j][i] = p;
        }
    }

    return;
}

void normalize(double** eigvec_matrix){

    /* Small routine to normalize a 4 4 matrix */

    int i, j;
    double norm = 0.0;

    for (i=0; i<4; i++){
        norm = 0.0;
        for (j=0; j<4; j++){
            norm += eigvec_matrix[i][j] * eigvec_matrix[i][j];
        }
        norm = sqrt(norm);
        for (j=0; j<4; j++){
            eigvec_matrix[i][j] /= norm;
        }
    }

    return;
}

void get_angle_quaternion(double* max_eigvec, double *angle){

    /* Get the tilt angle from the optimal rotation quaternion */

    /* Initialize stuff
       we won't use x and y */
    //~ double w = 0.0, x = 0.0, y = 0.0, z = 0.0;
    double w = 0.0, z = 0.0;

    /* Decompose the quaternion */
    w = max_eigvec[0];
    //~ x = max_eigvec[1];
    //~ y = max_eigvec[2];
    z = max_eigvec[3];

    /* Arctan2 to get the angle */
    /* The spin angle is 2 * (180/pi) * atan2(axis*q | q0)
       axis being [0, 0, 1] and q [x, y, z], axis*q = 1*z */
    *angle = (180.0/M_PI) * 2.0 * atan2(z, w);
    
    /* Get the angle between -180 and 180 */
    while (*angle < 180){
        *angle += 360;
    }
    while (*angle > 180){
        *angle -= 360;
    }

    return;

}

static PyMethodDef module_methods[] = {
    {"get_spinangle_traj",get_spinangle_traj,METH_VARARGS,get_spinangle_traj_docstring},
    {NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initspin(void)
{
    PyObject *m = Py_InitModule3("spin",module_methods,module_docstring);
    if (m == NULL)
        return;
    /* Load Numpy */
    import_array();
}
