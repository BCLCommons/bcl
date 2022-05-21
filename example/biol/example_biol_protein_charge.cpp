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
#include "biol/bcl_biol_protein_charge.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_protein_charge.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Jun 6, 2013
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolProteinCharge :
    public ExampleInterface
  {
  public:

    ExampleBiolProteinCharge *Clone() const
    {
      return new ExampleBiolProteinCharge( *this);
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
      io::IFStream read;

      // read seq1 from fasta
      BCL_Message( util::Message::e_Standard, "read fasta: " + AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from sequence
      biol::ProteinCharge protein_charge( seq1);

    /////////////////
    // data access //
    /////////////////

      // set the pk proerty to use
      protein_charge.SetPKProperty( biol::AATypeData::e_pK_Bjellqvist);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      const double pH( 7.0);
      storage::Map< biol::AATypeData::PropertyTypeEnum, double> charges;
      for
      (
        biol::AATypeData::PropertyTypeEnum itr( biol::ProteinCharge::s_FirstpKAAProperty);
        itr <= biol::ProteinCharge::s_LastpKAAProperty;
        ++itr
      )
      {
        BCL_MessageStd( "current property: " + itr.GetString());
        protein_charge.SetPKProperty( itr);
        charges[ itr] = protein_charge( pH);
      }

      util::GetLogger() << charges;

      BCL_ExampleCheckWithinAbsTolerance( charges[ biol::AATypeData::e_pK_Bjellqvist], 1.04509, 0.00009);
      BCL_ExampleCheckWithinAbsTolerance( charges[ biol::AATypeData::e_pK_ProMoST], 2.10073, 0.00009);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for biol::ProteinCharge");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( protein_charge);
      BCL_MessageVrb( "read object");
      biol::ProteinCharge protein_charge_read;
      ReadBCLObject( protein_charge_read);

      BCL_ExampleIndirectCheckWithinAbsTolerance( protein_charge_read( pH), protein_charge( pH), 0.00009, "read and write");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolProteinCharge

  const ExampleClass::EnumType ExampleBiolProteinCharge::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolProteinCharge())
  );

} // namespace bcl
