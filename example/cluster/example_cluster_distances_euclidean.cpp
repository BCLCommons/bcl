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
#include "cluster/bcl_cluster_distances_euclidean.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_distances_euclidean.cpp
  //!
  //! @author alexanns
  //! @date Jun 2, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterDistancesEuclidean :
    public ExampleInterface
  {
  public:

    ExampleClusterDistancesEuclidean *Clone() const
    {
      return new ExampleClusterDistancesEuclidean( *this);
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

      const linal::Vector< float> vector_a( 5, float( 1));
      const linal::Vector< float> vector_b( 5, float( 7));
      const float distance_correct( 13.4164);

      const util::SiPtr< const linal::Vector< float> > ptr_a( vector_a);
      const util::SiPtr< const linal::Vector< float> > ptr_b( vector_b);
      const storage::VectorND< 2, util::SiPtr< const linal::Vector< float> > > argument( ptr_a, ptr_b);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const cluster::DistancesEuclidean< float> calculator;

      // test clone
      util::ShPtr< cluster::DistancesEuclidean< float> > clone_constr( calculator.Clone());
      BCL_ExampleCheckWithinTolerance( clone_constr->operator()( argument), distance_correct, 0.001);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      const float distance( calculator( argument));

      BCL_MessageStd( "distance is " + util::Format()( distance));
      BCL_ExampleCheckWithinTolerance( distance, distance_correct, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

       // test Write and read
       WriteBCLObject( calculator);
       cluster::DistancesEuclidean< float> read_distance;
       ReadBCLObject( read_distance);
       BCL_ExampleCheckWithinTolerance( read_distance( argument), distance_correct, 0.001);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterDistancesEuclidean

  const ExampleClass::EnumType ExampleClusterDistancesEuclidean::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterDistancesEuclidean())
  );

} // namespace bcl
