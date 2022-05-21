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
#include "quality/bcl_quality_lcs.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_lcs.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityLCS :
    public ExampleInterface
  {
  public:

    ExampleQualityLCS *Clone() const
    {
      return new ExampleQualityLCS( *this);
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
      // get the protein model 1ubi
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel
        native( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone));

      // get the model for 1ubi
      const std::string model_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_rosetta_model.pdb"));
      assemble::ProteinModel
      model( Proteins::GetModel( model_filename, biol::GetAAClasses().e_AABackBone));

      // extract the CA atom coordinates from both native and model coordinates
      storage::Set< biol::AtomType> atom_types;
      atom_types.Insert( biol::GetAtomTypes().CA);
      util::SiPtrVector< const linal::Vector3D> native_coords( native.GetAtomCoordinates( atom_types));
      util::SiPtrVector< const linal::Vector3D> model_coords( model.GetAtomCoordinates( atom_types));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // initialize variables
      const double default_rmsd_cutoff( 5.0);
      const size_t default_seed_length( 3);
      const double rmsd_cutoff( 4.0);
      const size_t seed_length( 8);

      // test constructor with default values
      quality::LCS lcs_default;
      BCL_ExampleCheck( lcs_default.GetCutoff(), quality::LCS::GetDefaultRmsdCutoff());
      BCL_ExampleCheck( lcs_default.GetSeedLength(), quality::LCS::GetDefaultSeedLength());

      // construct with values
      quality::LCS lcs( rmsd_cutoff, seed_length);
      BCL_ExampleCheck( lcs.GetCutoff(), rmsd_cutoff);
      BCL_ExampleCheck( lcs.GetSeedLength(), seed_length);

      // clone construct
      util::ShPtr< quality::LCS> sp_lcs( lcs.Clone());
      BCL_ExampleCheck( sp_lcs->GetCutoff(), lcs.GetCutoff());
      BCL_ExampleCheck( sp_lcs->GetSeedLength(), lcs.GetSeedLength());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< quality::LCS>(), lcs.GetClassIdentifier());

      // check GetDefaultRmsdCutoff()
      BCL_ExampleCheck( quality::LCS::GetDefaultRmsdCutoff(), default_rmsd_cutoff);

      // check GetDefaultSeedLength()
      BCL_ExampleCheck( quality::LCS::GetDefaultSeedLength(), default_seed_length);

    ////////////////
    // operations //
    ////////////////

      BCL_MessageStd( "Testing CalculateMeasure()");

      // initialize expected values
      const size_t expected_lcs_1( 31);
      const size_t expected_lcs_2( 61);
      const size_t expected_lcs_5( 76);

      // calculate LCS for 1 angstrom cutoff
      const double lcs_1( quality::LCS( 1.0).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "LCS 1.0): " + util::Format()( lcs_1));
      BCL_ExampleCheck( lcs_1, expected_lcs_1);

      // calculate LCS for 2 angstrom cutoff
      const double lcs_2( quality::LCS( 2.0).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "LCS 2.0): " + util::Format()( lcs_2));
      BCL_ExampleCheck( lcs_2, expected_lcs_2);

      // calculate LCS for 5 angstrom cutoff
      const double lcs_5( quality::LCS( 5.0).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "LCS 5.0): " + util::Format()( lcs_5));
      BCL_ExampleCheck( lcs_5, expected_lcs_5);

      // expected superimposition matrix
      const double expected_superimposition_array[ 16] =
      {
         0.377536 , -0.910199,   0.170305   , 0,
        -0.0658895,  0.157044,   0.985391   , 0,
        -0.923648 , -0.383242,  -0.000682685, 0,
         0.210889 , 20.7056  , -31.0155     , 1
      };
      const linal::Vector< double> expected_superimposition_vector( 16, expected_superimposition_array);

      BCL_MessageStd( "testing CalculateSuperimposition");
      const linal::Matrix< double> superimposition_matrix
      (
        quality::LCS( 2.0).CalculateSuperimposition( native_coords, model_coords).GetMatrix()
      );
      BCL_MessageStd( "the superimposition matrix\n" + util::Format()( superimposition_matrix));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          linal::Vector< double>( 16, superimposition_matrix.Begin()), expected_superimposition_vector, 0.001, 0.001
        ), true, "comparison failed for transformation matrix calculated"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      WriteBCLObject( lcs);
      quality::LCS lcs_read;
      ReadBCLObject( lcs_read);
      BCL_ExampleCheck( lcs_read.GetCutoff(), lcs.GetCutoff());
      BCL_ExampleCheck( lcs_read.GetSeedLength(), lcs.GetSeedLength());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityLCS

  const ExampleClass::EnumType ExampleQualityLCS::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityLCS())
  );

} // namespace bcl
