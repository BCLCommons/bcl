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
#include "biol/bcl_biol_protein_params.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_protein_params.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date May 27, 2013
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolProteinParams :
    public ExampleInterface
  {
  public:

    ExampleBiolProteinParams *Clone() const
    {
      return new ExampleBiolProteinParams( *this);
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

      biol::ProteinParams protein_params;

      io::IFStream read;

      // read seq1 from fasta
      BCL_Message( util::Message::e_Standard, "read fasta: " + AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      BCL_Message( util::Message::e_Standard, "the fasta sequence is: ");
      seq1.WriteFasta( util::GetLogger());

      const storage::Table< double> params( protein_params( seq1));

      params.WriteFormatted( util::GetLogger());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolProteinParams

  const ExampleClass::EnumType ExampleBiolProteinParams::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolProteinParams())
  );

} // namespace bcl
