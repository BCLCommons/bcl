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
#include "biol/bcl_biol_ramachandran.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_ramachandran.cpp
  //!
  //! @author karakam
  //! @date Dec 30, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolRamachandran :
    public ExampleInterface
  {
  public:

    ExampleBiolRamachandran *Clone() const
    {
      return new ExampleBiolRamachandran( *this);
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

      // default constructor
      const biol::Ramachandran ramachandran( biol::Ramachandran::GetDefaultHistogramFilename());

      // create a reference on the static instance
      const biol::Ramachandran &default_ramachandran( biol::Ramachandran::GetDefaultInstance());

    /////////////////
    // data access //
    /////////////////

      // test GetDefaultSSTypeHistogramFilename()
      BCL_ExampleCheck( default_ramachandran.GetDefaultHistogramFilename(), "phi_psi_angles_by_sstype.histogram2D");

      // test GetSSTypeHistogramFilename()
      BCL_ExampleCheck( ramachandran.GetHistogramFilename(), "phi_psi_angles_by_sstype.histogram2D");

    ////////////////
    // operations //
    ////////////////

      // initialize range
      math::Range< double> range( -math::g_Pi, math::g_Pi);

      // initialize success to true
      bool success( true);

      BCL_MessageStd( "testing GetRandomPhiPsi for each AAType");
      // iterate over 20 aatypes
      for
      (
        biol::AATypes::const_iterator aatype_itr( biol::GetAATypes().Begin()),
          aatype_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
        aatype_itr != aatype_itr_end; ++aatype_itr
      )
      {
        // get a random number with this AAType
        const storage::VectorND< 2, double> phi_psi( ramachandran.GetRandomPhiPsi( *aatype_itr));

        // print the aatype and the random phi_psi
        BCL_MessageStd
        (
          ( *aatype_itr)->GetThreeLetterCode() + " ==> phi: " + util::Format()( phi_psi.First()) +
          "\tpsi: " + util::Format()( phi_psi.Second())
        );

        // make sure it is within valid range
        success &= range.IsWithin( phi_psi.First()) && range.IsWithin( phi_psi.Second());
      }

      // make sure all the phi psi values were within range
      BCL_ExampleIndirectCheck( success, true, "One or more random phi/psi generated from AAType were out of range");

      // set success back to true
      success = true;

      BCL_MessageStd( "testing GetRandomPhiPsi for each SSType and AAType");
      // iterate over SSTypes
      for
      (
        biol::SSTypes::const_iterator
        sstype_itr( biol::GetSSTypes().Begin()), sstype_itr_end( biol::GetSSTypes().COIL.GetIterator() + 1);
        sstype_itr != sstype_itr_end; ++sstype_itr
      )
      {
        // print the SSType
        BCL_MessageStd( "");
        // iterate over 20 aatypes
        for
        (
          biol::AATypes::const_iterator aatype_itr( biol::GetAATypes().Begin()),
            aatype_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
          aatype_itr != aatype_itr_end; ++aatype_itr
        )
        {
          // get a random number with this AAType
          const storage::VectorND< 2, double> phi_psi( ramachandran.GetRandomPhiPsi( *aatype_itr, *sstype_itr));

          // print the aatype and the random phi_psi
          BCL_MessageStd
          (
            sstype_itr->GetName() + "\t" + ( *aatype_itr)->GetThreeLetterCode() + " ==> phi: "
            + util::Format()( phi_psi.First()) + "\tpsi: " + util::Format()( phi_psi.Second())
          );

          // make sure it is within valid range
          success &= range.IsWithin( phi_psi.First()) && range.IsWithin( phi_psi.Second());
        }
      }

      // make sure all the phi psi values were within range
      BCL_ExampleIndirectCheck
      (
        success, true, "One or more random phi/psi generated from SSType and AATypes were out of range"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( ramachandran);

      // read back in
      biol::Ramachandran ramachandran_read;
      ReadBCLObject( ramachandran_read);

      // compare the histogram filename
      BCL_ExampleCheck( ramachandran.GetHistogramFilename(), ramachandran_read.GetHistogramFilename());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolRamachandran

  const ExampleClass::EnumType ExampleBiolRamachandran::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolRamachandran())
  );

} // namespace bcl
