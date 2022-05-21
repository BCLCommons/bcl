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
#include "contact/bcl_contact_recovery.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_recovery.cpp
  //!
  //! @author karakam
  //! @date Sep 22, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactRecovery :
    public ExampleInterface
  {
  public:

    ExampleContactRecovery *Clone() const
    {
      return new ExampleContactRecovery( *this);
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

      // initialize filenames
      const std::string filename_native( AddExampleInputPathToFilename( e_Biology, "2VQ4A_dssp.pdb"));
      const std::string filename_model_a( AddExampleInputPathToFilename( e_Biology, "2VQ4A_model1.pdb"));
      const std::string filename_model_b( AddExampleInputPathToFilename( e_Biology, "2VQ4A_model2.pdb"));

      // read the pdbs
      assemble::ProteinModel native_model( Proteins::GetModel( filename_native));
      assemble::ProteinModel model_a( Proteins::GetModel( filename_model_a));
      assemble::ProteinModel model_b( Proteins::GetModel( filename_model_b));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      contact::Recovery recovery_default;

      // constructor from a neighbor generator
      contact::Recovery recovery;

    ////////////////
    // operations //
    ////////////////

      // initialize expected contingency matrices
      const math::ContingencyMatrix expected_matrix_a( 134, 370, 348, 0);
      const math::ContingencyMatrix expected_matrix_b(  72, 242, 410, 0);

      // calculate the neighbor list containers for all models
      const assemble::AANeighborListContainer native_neigh
      (
        native_model.GetAminoAcids(), contact::g_ContactCbDistanceCutoff, contact::g_ContactMinSequenceSeparation, false
      );
      const assemble::AANeighborListContainer model_a_neigh
      (
        model_a.GetAminoAcids(), contact::g_ContactCbDistanceCutoff, contact::g_ContactMinSequenceSeparation, false
      );
      const assemble::AANeighborListContainer model_b_neigh
      (
        model_b.GetAminoAcids(), contact::g_ContactCbDistanceCutoff, contact::g_ContactMinSequenceSeparation, false
      );

      // calculate contingency matrix for native_model vs model_a neighbor list containers
      const math::ContingencyMatrix matrix_a( recovery.CalculateContingencyMatrix( native_neigh, model_a_neigh));
      BCL_MessageStd( "contingency_matrix_a:\n" + util::Format()( matrix_a));

      // calculate contingency matrix for native_model vs model_b neighbor list containers
      const math::ContingencyMatrix matrix_b( recovery.CalculateContingencyMatrix( native_neigh, model_b_neigh));
      BCL_MessageStd( "contingency_matrix_b:\n" + util::Format()( matrix_b));

    ///////////////
    // operators //
    ///////////////

      // initialize expected values
      const double expected_recovery_value_a( 100.0 * 134.0 / 482.0);
      const double expected_recovery_value_b( 100.0 *  72.0 / 482.0);

      // test operator with native_model itself
      BCL_MessageStd( "test operator() with native_model");
      const double recovery_value_native( recovery( native_model, native_model));
      BCL_MessageStd( "recovery_value_native: " + util::Format()( recovery_value_native));
      BCL_ExampleCheck( recovery_value_native, 100.0);

      // test operator with model a
      BCL_MessageStd( "test operator() with model_a");
      const double recovery_value_a( recovery( native_model, model_a));
      BCL_MessageStd( "recovery_value_a: " + util::Format()( recovery_value_a));
      BCL_ExampleCheck( recovery_value_a, expected_recovery_value_a);

      // test operator with model b
      BCL_MessageStd( "test operator() with model_b");
      const double recovery_value_b( recovery( native_model, model_b));
      BCL_MessageStd( "recovery_value_b: " + util::Format()( recovery_value_b));
      BCL_ExampleCheck( recovery_value_b, expected_recovery_value_b);

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( recovery);

      // read object
      contact::Recovery recovery_read;
      ReadBCLObject( recovery_read);

      // check that read object gives the same value
      const double recovery_read_value( recovery_read( native_model, model_a));
      BCL_ExampleCheck( recovery_read_value, expected_recovery_value_a);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactRecovery

  const ExampleClass::EnumType ExampleContactRecovery::s_Instance
  (
    GetExamples().AddEnum( ExampleContactRecovery())
  );
  
} // namespace bcl
