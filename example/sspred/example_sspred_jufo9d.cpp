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
#include "sspred/bcl_sspred_jufo9d.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sspred_jufo9d.cpp
  //!
  //! @author weinerbe
  //! @date May 14, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSspredJUFO9D :
    public ExampleInterface
  {
  public:

    ExampleSspredJUFO9D *Clone() const
    {
      return new ExampleSspredJUFO9D( *this);
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
      // attach read to a fasta file
      io::IFStream read;
      BCL_MessageStd( "read fasta: " + AddExampleInputPathToFilename( e_Biology, "2K73A.fasta"));
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "2K73A.fasta"));

      // create a sequence and read in the fasta sequence
      pdb::Factory factory( biol::GetAAClasses().e_AA);
      assemble::ProteinModel model( factory.ChainFromFastaStream( 'A', read));
      io::File::CloseClearFStream( read);

      // read the blast profile
      biol::BlastProfileHandler::TryReadProfileForProteinModel
      (
        model,
        AddExampleInputPathToFilename( e_Biology, "2K73A")
      );

    ////////////////
    // operations //
    ////////////////

      // calculate JUFO9D predictions for seq_a
      sspred::JUFO9D::Calculate( model);

    //////////////////////
    // input and output //
    //////////////////////

      // write out the predictions
      io::OFStream write;
      const std::string output_name( AddExampleOutputPathToFilename( sspred::GetNamespaceIdentifier(), "2K73A.jufo9d"));
      BCL_ExampleMustOpenOutputFile( write, output_name);
      sspred::MethodHandler::WritePredictionsForAASequence
      (
        write,
        *model.GetChains().FirstElement()->GetSequence(),
        sspred::GetMethods().e_JUFO9D,
        sspred::MethodHandler::e_NineState
      );
      io::File::CloseClearFStream( write);
      BCL_ExampleCheck( io::File::FilesMatch( output_name, output_name + ".correct"), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSspredJUFO9D

  const ExampleClass::EnumType ExampleSspredJUFO9D::s_Instance
  (
    GetExamples().AddEnum( ExampleSspredJUFO9D())
  );

} // namespace bcl
