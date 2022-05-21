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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "opencl/bcl_opencl_singular_value_decomposition.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_singular_value_decomposition.cpp
  //!
  //! @author loweew
  //! @date Aug 5, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclSingularValueDecomposition :
    public ExampleInterface
  {
  public:

    ExampleOpenclSingularValueDecomposition *Clone() const
    {
      return new ExampleOpenclSingularValueDecomposition( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

      // creating data set
      const size_t rows( 2);
      const size_t cols( 18);
      float data1[ 36] = {
                           1.0, 2.0, 3.0, 2.0, 2.0, 4.0, 3.0, 4.0, 3.0, 3.0, 4.0, 2.0, 3.0, 4.0, 3.0, 3.0, 4.0, 2.0,
                           4.0, 3.0, 2.0, 2.0, 3.0, 5.0, 5.0, 8.0, 7.0, 5.0, 2.0, 6.0, 4.0, 9.0, 1.0, 1.0, 4.0, 9.0
                         };
      linal::Matrix< float> sym_mat( rows, cols, data1);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construction
      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());
      opencl::SingularValueDecomposition< float> svd( queue);

      // running cpu standard
      linal::Matrix< float> eig_v( cols, cols), eig_u( rows, cols);
      linal::Matrix< float> eigenvalues( linal::SingularValueDecomposition( sym_mat, eig_v, eig_u));

      // running opencl svd
      svd( sym_mat);

      // checking results
      linal::Vector< double> eig_val_cpu( eigenvalues.GetNumberOfElements());
      linal::Vector< double> eig_val_gpu( eigenvalues.GetNumberOfElements());
      linal::Vector< double> eig_vec_u_cpu( eig_u.GetNumberOfElements());
      linal::Vector< double> eig_vec_u_gpu( eig_u.GetNumberOfElements());
      linal::Vector< double> eig_vec_v_cpu( eig_v.GetNumberOfElements());
      linal::Vector< double> eig_vec_v_gpu( eig_v.GetNumberOfElements());

      std::copy( eigenvalues.Begin(), eigenvalues.End(), eig_val_cpu.Begin());
      std::copy( svd.GetEigenValues().Begin(), svd.GetEigenValues().End(), eig_val_gpu.Begin());
      std::copy( eig_v.Begin(), eig_v.End(), eig_vec_v_cpu.Begin());
      std::copy( svd.GetEigenVectorsV().Begin(), svd.GetEigenVectorsV().End(), eig_vec_v_gpu.Begin());
      std::copy( eig_u.Begin(), eig_u.End(), eig_vec_u_cpu.Begin());
      std::copy( svd.GetEigenVectorsU().Begin(), svd.GetEigenVectorsU().End(), eig_vec_u_gpu.Begin());

      BCL_ExampleCheckWithinTolerance( eig_val_cpu  , eig_val_gpu  , 0.001);
      BCL_ExampleCheckWithinTolerance( eig_vec_v_cpu, eig_vec_v_gpu, 0.001);
//      BCL_ExampleCheckWithinTolerance( eig_vec_u_cpu, eig_vec_u_gpu, 0.001);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclSingularValueDecomposition

  const ExampleClass::EnumType ExampleOpenclSingularValueDecomposition::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclSingularValueDecomposition())
  );

} // namespace bcl
