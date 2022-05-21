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
#include "linal/bcl_linal_distance_geometry.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_distance_geometry.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author teixeipl
  //! @date Jun 8, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalDistanceGeometry :
    public ExampleInterface
  {
  public:

    ExampleLinalDistanceGeometry *Clone() const
    {
      return new ExampleLinalDistanceGeometry( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // Create a distance geometry function object
      linal::DistanceGeometry distance_geometry;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

// Data provided by David Nannemann's python distance geometry script found in ~nannemdp/scripts/distance_geometry.py
// applied to tetrahydrofuran molecule
//
//      original distance matrix
//      [[ 0.     1.428  1.531  2.405  2.285]
//       [ 1.428  0.     2.38   2.405  1.446]
//       [ 1.531  2.38   0.     1.528  2.429]
//       [ 2.405  2.405  1.528  0.     1.542]
//       [ 2.285  1.446  2.429  1.542  0.   ]]
//      eigenvalues
//      [  1.122e+01   3.999e+00   1.655e-01   0.000e+00  -8.970e-04]
//      eigenvectors
//      [[ 0.     0.     0.     1.     0.   ]
//       [ 0.236  0.574 -0.768  0.    -0.157]
//       [ 0.267 -0.614 -0.491  0.     0.558]
//       [ 0.688 -0.344  0.084  0.    -0.634]
//       [ 0.633  0.419  0.402  0.     0.512]]
//      output coordinates
//      [[ 0.     0.     0.     0.     0.   ]
//       [ 0.79   1.148 -0.313  0.    -0.   ]
//       [ 0.893 -1.227 -0.2    0.     0.   ]
//       [ 2.304 -0.688  0.034  0.    -0.   ]
//       [ 2.12   0.838  0.164  0.     0.   ]]
//      final_score 0.0
//      final distance matrix
//      [[ 0.     1.428  1.531  2.405  2.285]
//       [ 1.428  0.     2.38   2.405  1.446]
//       [ 1.531  2.38   0.     1.528  2.429]
//       [ 2.405  2.405  1.528  0.     1.542]
//       [ 2.285  1.446  2.429  1.542  0.   ]]

      // Create super simple test matrix
      storage::SymmetricMatrix< double> simple_distance_matrix( 5);

      // Set values within distance matrix
      simple_distance_matrix( 1, 0) = 1.428;
      simple_distance_matrix( 2, 0) = 1.531;
      simple_distance_matrix( 3, 0) = 2.405;
      simple_distance_matrix( 4, 0) = 2.285;
      simple_distance_matrix( 2, 1) = 2.38;
      simple_distance_matrix( 3, 1) = 2.405;
      simple_distance_matrix( 4, 1) = 1.446;
      simple_distance_matrix( 3, 2) = 1.528;
      simple_distance_matrix( 4, 2) = 2.429;
      simple_distance_matrix( 4, 3) = 1.542;

      // Debug for checking input is set up properly
      BCL_MessageDbg( "More precise distance matrix result:\n" + util::Format()( simple_distance_matrix));

      // Calculate the distance geometry for the more complex distance matrix specified above
      storage::Vector< linal::VectorND< double, 3> > simple_result( distance_geometry( simple_distance_matrix));

      // Print out matrix during debug
      BCL_MessageDbg( "Simple distance matrix result:\n" + util::Format()( simple_result));

      // Check output coordinates are right according to external calculation for super simple case (Nannemann data)
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 0)( 0)), 0.0,       0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 0)( 1)), 0.0,       0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 0)( 2)), 0.0,       0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 1)( 0)), 0.789937,  0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 1)( 1)), 1.14782,   0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 1)( 2)), 0.312579,  0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 2)( 0)), 0.893311,  0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 2)( 1)), 1.22734,   0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 2)( 2)), 0.199661,  0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 3)( 0)), 2.30426,   0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 3)( 1)), 0.68819,   0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 3)( 2)), 0.0339748, 0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 4)( 0)), 2.11969,   0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 4)( 1)), 0.837604,  0.001);
      BCL_ExampleCheckWithinTolerance( math::Absolute( simple_result( 4)( 2)), 0.163698,  0.001);

// TODO: change to indirect example check

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalDistanceGeometry

  const ExampleClass::EnumType ExampleLinalDistanceGeometry::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalDistanceGeometry())
  );

} // namespace bcl
