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
#include "biol/bcl_biol_blast_profile_handler.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_blast_profile_handler.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolBlastProfileHandler :
    public ExampleInterface
  {
  public:

    ExampleBiolBlastProfileHandler *Clone() const
    {
      return new ExampleBiolBlastProfileHandler( *this);
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
      // initialize streams
      io::IFStream read;
      io::OFStream write;

      // initialize filename
      const std::string
        fasta_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.fasta")),
        pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb")),
        ascii_filename( AddExampleInputPathToFilename( e_Biology, "1IE9A.ascii"));

      // read sequence from fasta
      BCL_MessageStd( "read fasta: " + fasta_filename);
      BCL_ExampleMustOpenInputFile( read, fasta_filename);
      biol::AASequence sequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read proteinmodel from pdb
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_MessageStd
      (
        "This class has the following identifier" + biol::BlastProfileHandler().GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // initialize a blast profile handler
      biol::BlastProfileHandler handler;

      // read blast profile for sequence
      BCL_MessageStd( "reading BLAST for sequence");
      BCL_ExampleMustOpenInputFile( read, ascii_filename);

      biol::BlastProfileHandler::ReadProfileForAASequence( read, sequence);
      io::File::CloseClearFStream( read);

      // tie ofstream to output file
      std::string output_filename( AddExampleOutputPathToFilename( handler, "generated_seq_1IE9.ascii"));
      BCL_ExampleMustOpenOutputFile( write, output_filename);

      // write profile for the sequence to a file
      BCL_MessageStd( "Writing Blast Profile for this sequence");
      biol::BlastProfileHandler::WriteProfileForAASequence( write, sequence);

      // read blast profile for the protein model
      BCL_MessageStd( "reading BLAST for protein model");
      BCL_ExampleAssert
      (
        biol::BlastProfileHandler::TryReadProfileForProteinModel
        (
          protein_model,
          AddExampleInputPathToFilename( e_Biology, "1IE9")
        ),
        true
      );

      // write profile for the sequence to a file
      BCL_MessageStd( "Writing Blast Profile for the protein model 1IE9");
      BCL_ExampleCheck
      (
        biol::BlastProfileHandler::WriteProfileForProteinModel
        (
          protein_model,
          "generated_model_1IE9",
          AddExampleOutputPathToFilename( handler, "")
        ),
        true
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolBlastProfileHandler

  const ExampleClass::EnumType ExampleBiolBlastProfileHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolBlastProfileHandler())
  );

} // namespace bcl
