static PyObject *get_spinangle_traj(PyObject *self, PyObject *args);
double get_spinangle(double *coords, double *ref, int num_atoms, double debug);
void build_matrix(double *coords, double *ref, double** quat_matrix, int num_atoms, double debug);
void jacobi(double** quaternion_matrix, double* eigval_matrix, double** eigvec_matrix, int num_atoms, int numb_rot);
void rotate_stuff(double** a, int i, int j, int k, int l, double* g, double* h, double tau, double s);
void eigsrt(double* eigval_matrix, double** eigvec_matrix);
void transpose(double** eigvec_matrix);
void normalize(double** eigvec_matrix);
void get_angle_quaternion(double* max_eigvec, double *angle);
