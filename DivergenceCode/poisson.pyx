import numpy
cimport numpy
ctypedef numpy.float_t DTYPE_f
ctypedef numpy.int_t DTYPE_i
cimport cython
@cython.boundscheck(False)

def poisson_iteration(numpy.ndarray[DTYPE_f, ndim=2] scalar,numpy.ndarray[DTYPE_i, ndim=2] coors,numpy.ndarray[DTYPE_f,ndim=2] div):
    #dims=scalar.shape
    cdef int dims[2]
    dims[0]=scalar.shape[0]
    dims[1]=scalar.shape[1]
    #cdef numpy.ndarray[int,ndim =1] dims = scalar.shape
    #scalar2 = numpy.zeros([dims[0],dims[1]])
    cdef numpy.ndarray scalar2 = numpy.zeros([dims[0], dims[1]])
    #cdef numpy.ndarray cdiv=div
    #cdef scalar2 = numpy.zeros_like(scalar)
    #print dims
    cdef int topj
    cdef int bottomj
    cdef int counter
    cdef int i
    cdef int j
    cdef double delta
    for counter in range(100):
        for i,j in coors:
            topj=(j+1+dims[1])%dims[1]
            bottomj=(j-1+dims[1])%dims[1]
            delta=0.25*(scalar[i,topj]+scalar[i,bottomj]+scalar[i+1,j]+scalar[i-1,j]-4*scalar[i,j]+4*div[i,j])
            scalar2[i,j]=scalar[i,j]+delta
            
        scalar=scalar2
        print counter
    return scalar

