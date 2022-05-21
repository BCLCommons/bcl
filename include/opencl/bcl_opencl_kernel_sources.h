// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_OPENCL_KERNEL_SOURCES_H_
#define BCL_OPENCL_KERNEL_SOURCES_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_kernel_source_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically
#include <map>
#include <CL/cl.hpp>

namespace bcl
{
  namespace opencl
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KernelSources
    //! @brief this class enumerates instance of KernelSourceInterface derived classes, that are used in different
    //!        algorithmic implementations that can be used in different opencl ports of general algorithms
    //!
    //! @see @link example_opencl_kernel_sources.cpp @endlink
    //! @author woetzen, loweew
    //! @date Jul 13, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KernelSources :
      public util::Enumerate< util::ShPtr< KernelSourceInterface>, KernelSources>
    {
      friend class util::Enumerate< util::ShPtr< KernelSourceInterface>, KernelSources>;
    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief parameter to change opencl kernels path
      static util::ShPtr< command::FlagInterface> &GetKernelsSourcePathFlag();

      //! @brief parameter to change opencl binary kernels path
      static util::ShPtr< command::FlagInterface> &GetKernelsBinaryPathFlag();

    //////////
    // data //
    //////////

      const KernelSource e_EuclideanDistance;             //!< creates pairwise distance matrix from two input matrices
      const KernelSource e_SvmKernelFunctions;            //!< svm kernel functions
      const KernelSource e_InsertionSort;                 //!< insertion sort kernel
      const KernelSource e_RMSD;                          //!< rmsd kernel
      const KernelSource e_SequentialMinimalOptimization; //!< smo kernel
      const KernelSource e_KappaNearestNeighbor;          //!< knn kernel
      const KernelSource e_ArtificialNeuralNetwork;       //!< ann kernel
      const KernelSource e_Linal;                         //!< linear algebra kernels
      const KernelSource e_Quality;                       //!< quality kernels
      const KernelSource e_Saxs;                          //!< SAXS kernels
      const KernelSource e_Buffer;                        //!< buffer kernels
      const KernelSource e_Development;                   //!< test kernels
      const KernelSource e_ClusterCenters;                //!<

    private:

      typedef std::map< std::string, std::map< CommandQueue, cl::Program> > ProgramMap;
      static ProgramMap &GetBuildPrograms();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      KernelSources();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add the Opencl kernel path to filename
      //! @param FILE_NAME filename of kernel
      //! @return kernel_path/FILE_NAME
      static std::string AddKernelPath( const std::string &FILE_NAME);

      //! @brief compile a kernel and return a binary
      //! @param KERNEL the kernel to compile
      //! @param PRECISION the precision to use
      //! @param QUEUE the command queue to program will run in
      //! @param OPTIONS compiler options
      //! @param ERROR_PTR location to store the ERROR at
      //! @return the cl::Program, either the stored already compiled version, or the newly compiled version
      static cl::Program
      Compile
      (
        const KernelSource &KERNEL,
        const util::CPPDataTypes::Types &PRECISION,
        const CommandQueue &QUEUE,
        const std::string &OPTIONS,
        cl_int *ERROR_PTR
      );

    }; // class KernelSources

    //! @brief construct on access function for all KernelSources
    //! @return reference to only instances of KernelSources
    BCL_API
    KernelSources &GetKernelSources();

  } // namespace opencl

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< opencl::KernelSourceInterface>, opencl::KernelSources>;

  } // namespace util
} // namespace bcl

#endif // BCL_OPENCL_KERNEL_SOURCES_H_
