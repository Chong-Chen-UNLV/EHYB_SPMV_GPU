/**
 * @file WildRiver.cpp
 * @brief Main function.
 * @author Dominique LaSalle <wildriver@domnet.org>
 * Copyright 2015-2016
 * @version 1
 */




#include <cstdint>
#include <iostream>
#include <memory>

#include "MatrixInHandle.hpp"
#include "MatrixOutHandle.hpp"
#include "GraphInHandle.hpp"
#include "GraphOutHandle.hpp"
#include "VectorInHandle.hpp"
#include "VectorOutHandle.hpp"
#include "Exception.hpp"




using namespace WildRiver;




/******************************************************************************
* HELPER FUNCTIONS ************************************************************
******************************************************************************/

namespace
{

/**
 * @brief Deleter function for using malloc with unique_ptr.
 */
struct c_delete
{
  void operator()(void * ptr)
  {
    free(ptr);
  }
};

}


/******************************************************************************
* API FUNCTIONS ***************************************************************
******************************************************************************/


extern "C" wildriver_matrix_handle * wildriver_open_matrix(
    char const * const filename,
    int const mode)
{
  try {
    std::unique_ptr<wildriver_matrix_handle> handle(
        new wildriver_matrix_handle);

    // initialize handle
    handle->mode = mode;
    handle->nrows = NULL_DIM;
    handle->ncols = NULL_DIM;
    handle->nnz = NULL_IND;
    handle->fd = nullptr;

    switch (mode) {
      case WILDRIVER_IN: {
          std::unique_ptr<MatrixInHandle> ptr(new MatrixInHandle(filename));
          ptr->getInfo(handle->nrows,handle->ncols,handle->nnz);
          handle->fd = reinterpret_cast<void*>(ptr.release());
        }
        break;
      case WILDRIVER_OUT: {
          std::unique_ptr<MatrixOutHandle> ptr(new MatrixOutHandle(filename));
          handle->fd = reinterpret_cast<void*>(ptr.release());
        }
        break;
      default:
        throw BadParameterException(std::string("Unknown matrix mode: ") + \
            std::to_string(mode));
    }

    return handle.release();
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to open matrix due to: " << e.what() \
        << std::endl;
    return nullptr;
  }
}


extern "C" int wildriver_load_matrix(
    wildriver_matrix_handle * const handle,
    ind_t * const rowptr,
    dim_t * const rowind,
    val_t * const rowval,
    double * progress)
{
  if (progress != nullptr) {
    *progress = 0.0;
  }

  try {
    if (handle->mode != WILDRIVER_IN) {
      throw BadParameterException( \
          std::string("Cannot load matrix in mode: ") + \
          std::to_string(handle->mode));
    }

    // check input
    if (handle->nrows == NULL_DIM) {
      throw BadParameterException("Number of rows has not been set.");
    }
    if (handle->ncols == NULL_DIM) {
      throw BadParameterException("Number of columns has not been set.");
    }
    if (handle->nnz == NULL_IND) {
      throw BadParameterException("Number of non-zeros has not been set.");
    }
    if (handle->fd == nullptr) {
      throw BadParameterException("The file descriptor has not been set.");
    }

    MatrixInHandle * const inHandle = \
        reinterpret_cast<MatrixInHandle*>(handle->fd);

    // allocate matrix
    inHandle->readSparse(rowptr,rowind,rowval,progress);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to read matrix due to: " << e.what() \
        << std::endl;
    return 0;
  }

  if (progress != nullptr) {
    *progress = 1.0;
  }

  return 1;
}


extern "C" int wildriver_save_matrix(
    wildriver_matrix_handle * const handle,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval,
    double * const progress)
{
  if (progress != nullptr) {
    *progress = 0.0;
  }

  try {
    if (handle->mode != WILDRIVER_OUT) {
      throw BadParameterException( \
          std::string("Cannot save matrix in mode: ") + \
          std::to_string(handle->mode));
    }

    // check input
    if (handle->nrows == NULL_DIM) {
      throw BadParameterException("Number of rows has not been set.");
    }
    if (handle->ncols == NULL_DIM) {
      throw BadParameterException("Number of columns has not been set.");
    }
    if (handle->nnz == NULL_IND) {
      throw BadParameterException("Number of non-zeros has not been set.");
    }
    if (handle->fd == nullptr) {
      throw BadParameterException("The file descriptor has not been set.");
    }

    MatrixOutHandle * const outHandle = \
        reinterpret_cast<MatrixOutHandle*>(handle->fd);

    // save matrix
    outHandle->setInfo(handle->nrows,handle->ncols,handle->nnz);
    outHandle->writeSparse(rowptr,rowind,rowval);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to save matrix due to: " << e.what() \
        << std::endl;
    return 0;
  }

  if (progress != nullptr) {
    *progress = 1.0;
  }

  return 1;
}


extern "C" void wildriver_close_matrix(
    wildriver_matrix_handle * handle)
{
  if (handle->fd != nullptr) {
    if (handle->mode == WILDRIVER_IN) {
      delete reinterpret_cast<MatrixInHandle*>(handle->fd);
    } else if (handle->mode == WILDRIVER_OUT) {
      delete reinterpret_cast<MatrixOutHandle*>(handle->fd);
    } else {
      std::cerr << "WARNING: unable to determine current mode" << \
          handle->mode << std::endl;
    }
  }

  delete handle;
}


extern "C" wildriver_vector_handle * wildriver_open_vector(
    char const * const filename,
    int const mode)
{
  try {
    std::unique_ptr<wildriver_vector_handle> handle(
        new wildriver_vector_handle);

    // initialize handle
    handle->mode = mode;
    handle->size = NULL_IND;
    handle->fd = nullptr;

    switch (mode) {
      case WILDRIVER_IN: {
          std::unique_ptr<VectorInHandle> ptr(new VectorInHandle(filename));
          handle->size = ptr->getSize();
          handle->fd = reinterpret_cast<void*>(ptr.release());
        }
        break;
      case WILDRIVER_OUT: {
          std::unique_ptr<VectorOutHandle> ptr(new VectorOutHandle(filename));
          handle->fd = reinterpret_cast<void*>(ptr.release());
        }
        break;
      default:
        throw BadParameterException(std::string("Unknown vector mode: ") + \
            std::to_string(mode));
    }

    return handle.release();
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to open vector due to: " << e.what() \
        << std::endl;
    return nullptr;
  }
}


extern "C" int wildriver_load_vector(
    wildriver_vector_handle * handle,
    val_t * const vals,
    double * progress)
{
  if (progress != nullptr) {
    *progress = 0.0;
  }

  try {
    if (handle->mode != WILDRIVER_IN) {
      throw BadParameterException( \
          std::string("Cannot load vector in mode: ") + \
          std::to_string(handle->mode));
    }

    // check input
    if (handle->size == NULL_IND) {
      throw BadParameterException("The size of the vector has not been set.");
    }
    if (handle->fd == nullptr) {
      throw BadParameterException("The file descriptor has not been set.");
    }

    VectorInHandle * const inHandle = \
        reinterpret_cast<VectorInHandle*>(handle->fd);

    // allocate vector
    inHandle->read(vals,progress);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to read vector due to: " << e.what() \
        << std::endl;
    return 0;
  }

  if (progress != nullptr) {
    *progress = 1.0;
  }

  return 1;
}


extern "C" int wildriver_save_vector(
    wildriver_vector_handle * handle,
    val_t const * const vals,
    double * progress)
{
  if (progress != nullptr) {
    *progress = 0.0;
  }

  try {
    if (handle->mode != WILDRIVER_OUT) {
      throw BadParameterException( \
          std::string("Cannot save vector in mode: ") + \
          std::to_string(handle->mode));
    }

    // check input
    if (handle->size == NULL_IND) {
      throw BadParameterException("Size of vector has not been set.");
    }
    if (handle->fd == nullptr) {
      throw BadParameterException("The file descriptor has not been set.");
    }

    VectorOutHandle * const outHandle = \
        reinterpret_cast<VectorOutHandle*>(handle->fd);

    // save matrix
    outHandle->write(vals,handle->size,progress);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to save vector due to: " << e.what() \
        << std::endl;
    return 0;
  }

  if (progress != nullptr) {
    *progress = 1.0;
  }

  return 1;
}


extern "C" void wildriver_close_vector(
    wildriver_vector_handle * handle)
{
  if (handle->fd != nullptr) {
    if (handle->mode == WILDRIVER_IN) {
      delete reinterpret_cast<VectorInHandle*>(handle->fd);
    } else if (handle->mode == WILDRIVER_OUT) {
      delete reinterpret_cast<VectorOutHandle*>(handle->fd);
    } else {
      std::cerr << "WARNING: unable to determine current mode" << \
          handle->mode << std::endl;
    }
  }

  delete handle;
}




/******************************************************************************
* DEPRECATED FUNCTIONS ********************************************************
******************************************************************************/


extern "C" int wildriver_read_matrix(
    char const * const fname,
    dim_t * const r_nrows,
    dim_t * const r_ncols,
    ind_t * const r_nnz,
    ind_t ** const r_rowptr,
    dim_t ** const r_rowind,
    val_t ** const r_rowval)
{
  try {
    MatrixInHandle handle(fname);

    dim_t nrows, ncols;
    ind_t nnz;

    // read the header
    handle.getInfo(nrows,ncols,nnz);

    // allocate matrix
    size_t nbytes = sizeof(ind_t)*(nrows+1);
    std::unique_ptr<ind_t,c_delete> rowptr((ind_t*)malloc(nbytes));
    if (!rowptr.get()) {
      throw OutOfMemoryException(nbytes);
    }

    nbytes = sizeof(dim_t)*nnz;
    std::unique_ptr<dim_t,c_delete> rowind((dim_t*)malloc(nbytes));
    if (!rowind.get()) {
      throw OutOfMemoryException(nbytes);
    }

    std::unique_ptr<val_t,c_delete> rowval;
    if (r_rowval) {
      // we need to use rowval
      nbytes = sizeof(val_t)*nnz;
      rowval.reset((val_t*)malloc(nbytes));
      if (!rowind.get()) {
        throw OutOfMemoryException(nbytes);
      }
    } else {
      // don't allocate rowval
      nbytes = 0;
    }

    handle.readSparse(rowptr.get(),rowind.get(),rowval.get());

    // we've completely succeed -- assign pointers
    *r_nrows = nrows;
    *r_ncols = ncols;
    *r_nnz = nnz;

    *r_rowptr = rowptr.release();
    *r_rowind = rowind.release();
    if (r_rowval) {
      *r_rowval = rowval.release();
    }
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to read matrix due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}


extern "C" int wildriver_write_matrix(
    char const * const fname,
    dim_t const nrows,
    dim_t const ncols,
    ind_t const nnz,
    ind_t const * const rowptr,
    dim_t const * const rowind,
    val_t const * const rowval)
{
  try{
    MatrixOutHandle handle(fname);

    handle.setInfo(nrows,ncols,nnz);

    handle.writeSparse(rowptr,rowind,rowval);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to write matrix due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}


extern "C" int wildriver_read_graph(
    char const * const fname,
    dim_t * const r_nvtxs,
    ind_t * const r_nedges,
    int * const r_nvwgts,
    int * const r_ewgts,
    ind_t ** const r_xadj,
    dim_t ** const r_adjncy,
    val_t ** const r_vwgt,
    val_t ** const r_adjwgt)
{
  try {
    GraphInHandle handle(fname);

    dim_t nvtxs;
    ind_t nedges;
    int nvwgts;
    bool ewgts;

    // read the header
    handle.getInfo(nvtxs,nedges,nvwgts,ewgts);

    // allocate matrix
    size_t nbytes = sizeof(ind_t)*(nvtxs+1);
    std::unique_ptr<ind_t,c_delete> xadj((ind_t*)malloc(nbytes));
    if (!xadj.get()) {
      throw OutOfMemoryException(nbytes);
    }

    nbytes = sizeof(dim_t)*nedges;
    std::unique_ptr<dim_t,c_delete> adjncy((dim_t*)malloc(nbytes));
    if (!adjncy.get()) {
      throw OutOfMemoryException(nbytes);
    }

    std::unique_ptr<val_t> vwgt;
    if (r_vwgt && nvwgts > 0) {
      // we need to use rowval
      nbytes = sizeof(val_t)*nvtxs*nvwgts;
      vwgt.reset((val_t*)malloc(nbytes));
      if (!adjncy.get()) {
        throw OutOfMemoryException(nbytes);
      }
    }

    std::unique_ptr<val_t,c_delete> adjwgt;
    if (r_adjwgt) {
      // we need to use rowval
      nbytes = sizeof(val_t)*nedges;
      adjwgt.reset((val_t*)malloc(nbytes));
      if (!adjncy.get()) {
        throw OutOfMemoryException(nbytes);
      }
    }

    handle.readGraph(xadj.get(),adjncy.get(),vwgt.get(),adjwgt.get());

    // we've completed exception possible tasks -- assign pointers
    *r_xadj = xadj.release();
    *r_adjncy = adjncy.release();
    if (r_vwgt) {
      *r_vwgt = vwgt.release();
    }
    if (r_adjwgt) {
      *r_adjwgt = adjwgt.release();
    }

    *r_nvtxs = nvtxs;

    if (r_nedges) {
      *r_nedges = nedges;
    }

    if (r_ewgts) {
      // convert to c style int
      *r_ewgts = (int)ewgts;
    }
    if (r_nvwgts) {
      *r_nvwgts = nvwgts;
    }
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to read graph due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}


extern "C" int wildriver_write_graph(
    char const * const fname,
    dim_t const nvtxs,
    ind_t const nedges,
    int nvwgts,
    ind_t const * const xadj,
    dim_t const * const adjncy,
    val_t const * const vwgt,
    val_t const * const adjwgt)
{
  try {
    GraphOutHandle handle(fname);

    bool const ewgts = adjwgt != NULL;

    handle.setInfo(nvtxs,nedges,nvwgts,ewgts);

    handle.writeGraph(xadj,adjncy,vwgt,adjwgt);
  } catch (std::exception const & e) {
    std::cerr << "ERROR: failed to write graph due to: " << e.what() \
        << std::endl;
    return 0;
  }

  return 1;
}
