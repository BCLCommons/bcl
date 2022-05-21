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
#include "restraint/bcl_restraint_distance.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_distance.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintDistance :
    public ExampleInterface
  {
  public:

    ExampleRestraintDistance *Clone() const
    {
      return new ExampleRestraintDistance( *this);
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

      // initialize variable
      const double distance_value( 5.0);
      const double upper_bound( 8.0);
      const double lower_bound( 3.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::Distance distance_def;
      BCL_ExampleCheck( util::IsDefined( distance_def.GetDistance()), false);
      BCL_ExampleCheck( util::IsDefined( distance_def.UpperBound()), false);
      BCL_ExampleCheck( util::IsDefined( distance_def.LowerBound()), false);

      // test constructor from a distance and bounds
      restraint::Distance distance( distance_value, upper_bound, lower_bound);

      // test copy constructor
      restraint::Distance distance_copy( distance);
      BCL_ExampleCheck( distance_copy.GetDistance(), distance.GetDistance());
      BCL_ExampleCheck( distance_copy.UpperBound(), distance.UpperBound());
      BCL_ExampleCheck( distance_copy.LowerBound(), distance.LowerBound());

      // test clone
      util::ShPtr< restraint::Distance> sp_distance( distance.Clone());
      BCL_ExampleCheck( sp_distance->GetDistance(), distance.GetDistance());
      BCL_ExampleCheck( sp_distance->UpperBound(), distance.UpperBound());
      BCL_ExampleCheck( sp_distance->LowerBound(), distance.LowerBound());

    /////////////////
    // data access //
    /////////////////

      // test the static class name
      BCL_ExampleCheck( GetStaticClassName( distance), "bcl::restraint::Distance");

      // test GetClassIdentifier()
      BCL_ExampleCheck( distance.GetClassIdentifier(), "bcl::restraint::Distance");

      // test GetIdentification()
      BCL_ExampleCheck( distance.GetIdentification(), "3 5 8");

      // test GetDistance()
      BCL_ExampleCheck( distance.GetDistance(), distance_value);

      // test UpperBound()
      BCL_ExampleCheck( distance.UpperBound(), upper_bound);

      // test LowerBound()
      BCL_ExampleCheck( distance.LowerBound(), lower_bound);

      // test IsDefined()
      BCL_ExampleCheck( distance.IsDefined(), true);
      BCL_ExampleCheck( restraint::Distance().IsDefined(), false);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test write
      WriteBCLObject( distance);

      // create new object and test read
      restraint::Distance distance_read;
      ReadBCLObject( distance_read);

      BCL_ExampleIndirectCheck
      (
        distance.GetDistance() == distance_read.GetDistance() &&
        distance.UpperBound() == distance_read.UpperBound() &&
        distance.LowerBound() == distance_read.LowerBound(),
        true,
        "comparison for the written and read objects"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // test GetUpperError()
      BCL_ExampleCheck( distance.GetUpperError(), upper_bound - distance_value);

      // test GetLowerError()
      BCL_ExampleCheck( distance.GetLowerError(), distance_value - lower_bound);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintDistance

  const ExampleClass::EnumType ExampleRestraintDistance::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintDistance())
  );

} // namespace bcl
