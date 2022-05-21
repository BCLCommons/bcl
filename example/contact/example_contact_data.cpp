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
#include "contact/bcl_contact_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_data.cpp
  //!
  //! @author heinzes1
  //! @date Aug 5, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactData :
    public ExampleInterface
  {
  public:

    ExampleContactData *Clone() const
    {
      return new ExampleContactData( *this);
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

      // test default constructor
      contact::Data default_data;

      // test Clone below

    /////////////////
    // data access //
    /////////////////

      // test reading from MergedPrediction
      BCL_ExampleIndirectCheck( util::IsDefined( default_data.MergedPrediction()), false, "Read MergedPrediction");

      // test writing to MergedPrediction
      default_data.MergedPrediction() = math::g_Pi;
      BCL_ExampleIndirectCheck( default_data.MergedPrediction(), math::g_Pi, "Write MergedPrediction");

      // test IsDefined
      BCL_ExampleIndirectCheck( default_data.IsDefined(), false, "IsDefined is false");

    ////////////////
    // operations //
    ////////////////

      // test Swap below

    ///////////////
    // operators //
    ///////////////

      // test operator[]
      BCL_ExampleIndirectCheck
      (
        util::IsDefined( default_data[ contact::GetTypes().HELIX_STRAND]),
        false,
        "Read HELIX_STRAND prediction"
      );

      const double helix_strand_prediction( 1.0), strand_helix_prediction( 2.0);

      // test writing contact prediction values
      default_data[ contact::GetTypes().HELIX_STRAND] = helix_strand_prediction;
      default_data[ contact::GetTypes().STRAND_HELIX] = strand_helix_prediction;
      BCL_ExampleIndirectCheck
      (
        default_data[ contact::GetTypes().HELIX_STRAND],
        helix_strand_prediction,
        "Read HELIX_STRAND prediction"
      );

      // test Swap
      default_data.Swap( contact::GetTypes().HELIX_STRAND, contact::GetTypes().STRAND_HELIX);
      BCL_ExampleIndirectCheck
      (
        default_data[ contact::GetTypes().HELIX_STRAND],
        strand_helix_prediction,
        "Swap contact types"
      );

      // test Clone
      util::ShPtr< contact::Data> cloned_data( default_data.Clone());
      BCL_ExampleIndirectCheck( cloned_data->MergedPrediction(), math::g_Pi, "Cloned MergedPrediction");
      BCL_ExampleIndirectCheck
      (
        ( *cloned_data)[ contact::GetTypes().HELIX_STRAND],
        strand_helix_prediction,
        "Test clone read HELIX_STRAND"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // initialize all contact type in data, so that we do not compare nan
      double prediction( 3.0);
      for( size_t type_pos( 0); type_pos < contact::Types::s_NumberValidTypes; ++type_pos, ++prediction)
      {
        contact::Type type( contact::GetTypes().GetEnumFromIndex( type_pos));
        ( *cloned_data)[ type] = prediction;
      }

      // test Write and Read
      WriteBCLObject( *cloned_data);
      contact::Data empty_data;
      ReadBCLObject( empty_data);
      BCL_ExampleIndirectAssert( IsDataEqual( *cloned_data, empty_data), true, "Read and Write");

      // test IsDefined
      BCL_ExampleIndirectCheck( cloned_data->IsDefined(), true, "IsDefined is true");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    //! @brief test if two contact::Data are equal
    //! @param DATA_A first data to compare
    //! @param DATA_B second data to compare
    //! @return if DATA_A and DATA_B are equal
    bool IsDataEqual( const contact::Data &DATA_A, const contact::Data &DATA_B) const
    {
      // check if merged contact probability is equal
      if( !math::EqualWithinTolerance( DATA_A.MergedPrediction(), DATA_B.MergedPrediction()))
      {
        return false;
      }

      // check if prediction for all contact types are equal
      for( size_t type_pos( 0); type_pos < contact::Types::s_NumberValidTypes; ++type_pos)
      {
        contact::Type type( contact::GetTypes().GetEnumFromIndex( type_pos));
        if( !math::EqualWithinTolerance( DATA_A[ type], DATA_B[ type]))
        {
          return false;
        }
      }

      return true;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactData

  const ExampleClass::EnumType ExampleContactData::s_Instance( GetExamples().AddEnum( ExampleContactData()));
  
} // namespace bcl
