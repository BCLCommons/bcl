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
#include "release/bcl_app_jufo.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_jufo.cpp
  //!
  //! @author weinerbe
  //! @date September 17, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppJufo :
    public ExampleInterface
  {
  public:

    ExampleAppJufo *Clone() const
    {
      return new ExampleAppJufo( *this);
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

      ApplicationExampleHelper jufo_helper( app::Jufo::Jufo_Instance);

    ////////////////
    // operations //
    ////////////////

      // check that flags are needed
      BCL_ExampleCheck( jufo_helper.CheckCommandString( false), false);

      // get input and output file names
      const std::string input_file( AddExampleInputPathToFilename( e_Biology, "2K73A.fasta"));
      const std::string output_filename
      (
        AddExampleOutputPathToFilename( sspred::GetNamespaceIdentifier(), "1ubiA_app_test.jufo")
      );

      // add flags
      jufo_helper.ResetFlagsAndParameters();
      jufo_helper.AddParameter( input_file);

      // check the command line
      BCL_ExampleAssert( jufo_helper.CheckCommandString( true), true);

      // run the command line
      if( BCL_ExampleCheck( jufo_helper.RunCommand(), 0))
      {
        // read in the pdb
        assemble::ProteinModel protein_model
        (
          Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
        );

        // set membrane
        util::ShPtr< assemble::ProteinModelData> sp_data( protein_model.GetProteinModelData());
        sp_data->Insert
        (
          assemble::ProteinModelData::e_Membrane,
          util::ShPtr< biol::Membrane>( new biol::Membrane())
        );
        protein_model.SetProteinModelData( sp_data);

        // update PDB "prediction" w/ environment types
        sspred::PDB::SetEnvironmentTypes( protein_model);

        // get the Chain A sequence
        biol::AASequence seq( *( protein_model.GetChain( 'A')->GetSequence()));

        // read in the JUFO prediction
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "2K73A.jufo9d"));
        sspred::MethodHandler::ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_JUFO9D);
        io::File::CloseClearFStream( read);

        // initialize correct prediction counter
        size_t correct_predictions( 0);

        // iterate over the sequence
        for
        (
          util::ShPtrVector< biol::AABase>::const_iterator aa_itr( seq.Begin()), aa_itr_end( seq.End());
          aa_itr != aa_itr_end; ++aa_itr
        )
        {
          // if the JUFO prediction is the same as the PDB file
          if
          (
            ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_JUFO9D)->GetOneStateSSPrediction() ==
              ( *aa_itr)->GetSSPrediction( sspred::GetMethods().e_PDB)->GetOneStateSSPrediction()
          )
          {
            // increment the correct counter
            ++correct_predictions;
          }
        }

        // expect JUFO to correctly predict 142
        const size_t nr_expected_correct( 142);
        BCL_ExampleCheck( correct_predictions >= nr_expected_correct, true);
        BCL_MessageStd
        (
          "JUFO predicted " + util::Format()( correct_predictions) + " out of " +
            util::Format()( seq.GetSize()) + " correctly"
        );
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppJufo

  const ExampleClass::EnumType ExampleAppJufo::s_Instance
  (
    GetExamples().AddEnum( ExampleAppJufo())
  );

} // namespace bcl
