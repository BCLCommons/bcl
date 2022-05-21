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
#include "quality/bcl_quality_dmf.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "linal/bcl_linal_matrix.h"
#include "quality/bcl_quality_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_dmf.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityDMF :
    public ExampleInterface
  {
  public:

    ExampleQualityDMF *Clone() const
    {
      return new ExampleQualityDMF( *this);
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
      const storage::Set< double> cutoff_set( storage::Set< double>::Create( 2.0, 3.0, 4.0, 5.0));

      // default constructor
      quality::DMF dmf( quality::Measures::GetDistanceCutoffsTS());
      BCL_ExampleCheck( dmf.GetDistanceCutoffs(), quality::Measures::GetDistanceCutoffsTS());

      // constructor from vector
      quality::DMF dmf_b( cutoff_set);
      BCL_ExampleCheck( dmf_b.GetDistanceCutoffs(), cutoff_set);

      // clone
      util::ShPtr< quality::DMF> sp_dmf( dmf.Clone());
      BCL_ExampleCheck( sp_dmf->GetDistanceCutoffs(), dmf.GetDistanceCutoffs());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_ExampleCheck( GetStaticClassName< quality::DMF>(), dmf.GetClassIdentifier());

    ////////////////
    // operations //
    ////////////////

      // if the current message level is verbose
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
      {
        // iterate over different distances
        for( double cutoff( 0.5), end_cutoff( 10.0); cutoff <= end_cutoff; cutoff += 0.50)
        {
          // calculate the corresponding gdt with this cutoff
          const double this_dmf( quality::DMF::CalculateMeasure( native_coords, model_coords, cutoff));
          // print out the result
          BCL_MessageVrb( "DMF " + util::Format()( cutoff) + " : " + util::Format()( this_dmf));
        }
      }

      // initialize expected values
      const double expected_dmf_ts( 91.2544);
      const double expected_dmf_ha( 84.8596);
      const double expected_dmf_3( 92.6667);
      const double expected_dmf_matrix( 88.8889);
      const linal::Matrix< double> matrix
      (
        coord::CalculateDifferenceDistanceMatrix
        (
          native_coords.SubSiPtrVector( 0, 10), model_coords.SubSiPtrVector( 0, 10)
        )
      );

      // calculate DMF_TS which is ( DMF_1 + DMF_2 + DMF_4 + DMF_8) / 4
      BCL_MessageStd( "test operator() for DMF_TS");
      const double dmf_ts( ( *quality::GetMeasures().e_DMF_TS)->CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "DMF_TS: " + util::Format()( dmf_ts));

      // calculate DMF_HA which is ( DMF_0.5 + DMF_1 + DMF_2 + DMF_4) / 4
      BCL_MessageStd( "test operator() for DMF_HA");
      const double dmf_ha( ( *quality::GetMeasures().e_DMF_HA)->CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "DMF_HA: " + util::Format()( dmf_ha));

      // test operator with a cutoff of 3
      BCL_MessageStd( "test operator() with a cutoff of 3");
      const double dmf_3( dmf.CalculateMeasure( native_coords, model_coords, 3));
      BCL_MessageStd( "DMF_3: " + util::Format()( dmf_3));

      // test operator with a matrix
      BCL_MessageStd( "test operator() with a matrix");
      const double dmf_matrix( dmf.CalculateMeasureMatrixCutoff( matrix, 0.5));
      BCL_MessageStd( "DMF_matrix: " + util::Format()( dmf_matrix));

      // check the results
      BCL_ExampleCheckWithinTolerance( dmf_ts, expected_dmf_ts, 0.001);
      BCL_ExampleCheckWithinTolerance( dmf_ha, expected_dmf_ha, 0.001);
      BCL_ExampleCheckWithinTolerance( dmf_3, expected_dmf_3, 0.001);
      BCL_ExampleCheckWithinTolerance( dmf_matrix, expected_dmf_matrix, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      WriteBCLObject( dmf);
      quality::DMF dmf_read( storage::Set< double>( 0.0));
      ReadBCLObject( dmf_read);
      BCL_ExampleCheck( dmf_read.GetDistanceCutoffs(), dmf.GetDistanceCutoffs());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityDMF

  const ExampleClass::EnumType ExampleQualityDMF::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityDMF())
  );

} // namespace bcl
