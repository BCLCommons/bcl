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
#include "restraint/bcl_restraint_cone_model.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_ensemble.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_cone_model.cpp
  //!
  //! @author alexanns
  //! @date Mar 4, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintConeModel :
    public ExampleInterface
  {
  public:

    ExampleRestraintConeModel *Clone() const
    {
      return new ExampleRestraintConeModel( *this);
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
      // create proteins
      const util::ShPtr< assemble::ProteinModel> model_a
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "2LZM_sl_105_0001.pdb"), biol::GetAAClasses().e_AAComplete
        ).Clone()
      );
      const util::ShPtr< assemble::ProteinModel> model_b
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "2LZM_sl_105_0002.pdb"), biol::GetAAClasses().e_AAComplete
        ).Clone()
      );
      const util::ShPtr< assemble::ProteinModel> model_c
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "2LZM_sl_105_0003.pdb"), biol::GetAAClasses().e_AAComplete
        ).Clone()
      );

      // make and fill ensemble with proteins
      assemble::ProteinEnsemble ensemble;
      ensemble.InsertElement( model_a);
      ensemble.InsertElement( model_b);
      ensemble.InsertElement( model_c);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      restraint::ConeModel def_constr;

      // clone constructor
      util::ShPtr< restraint::ConeModel> clone_constr( def_constr.Clone());

    /////////////////
    // data access //
    /////////////////

      // GetStaticClassName
      const std::string correct_static_class_name( "bcl::restraint::ConeModel");
      BCL_ExampleCheck( GetStaticClassName< restraint::ConeModel>(), correct_static_class_name);

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< restraint::ConeModel>(), clone_constr->GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      { // SLCBSLMaxAngle
        const double result( restraint::ConeModel::SLCBSLMaxAngle( ensemble));
        BCL_MessageDbg( "calculated result " + util::Format()( result));
        const double expected( 0.841783);
        BCL_MessageDbg( "expected result " + util::Format()( expected));
        BCL_ExampleCheckWithinTolerance( result, expected, 0.000001);
      }

      { // SLeffectiveCBCAAngle
        const double result( restraint::ConeModel::SLeffectiveCBCAAngle( ensemble));
        BCL_MessageDbg( "calculated result " + util::Format()( result));
        const double expected( 1.82451);
        BCL_MessageDbg( "expected result " + util::Format()( expected));
        BCL_ExampleCheckWithinTolerance( result, expected, 0.000001);
      }

      { // SLeffectiveCBDistance
        const double result( restraint::ConeModel::SLeffectiveCBDistance( ensemble));
        BCL_MessageDbg( "calculated result " + util::Format()( result));
        const double expected( 6.49789);
        BCL_MessageDbg( "expected result " + util::Format()( expected));
        BCL_ExampleCheckWithinTolerance( result, expected, 0.000001);
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( def_constr);

      // read the object back in
      restraint::ConeModel read;
      ReadBCLObject( read);

    //////////////////////
    // helper functions //
    //////////////////////

      // GetSpinLabelResidue
      const util::SiPtr< const biol::AABase> sl( restraint::ConeModel::GetSpinLabelResidue( *model_a));
      BCL_ExampleCheck( sl->GetChainID(), 'A');
      BCL_ExampleCheck( sl->GetSeqID(), 105);

      { // GetUnpairedElectronCoordinates
        const linal::Vector3D result( restraint::ConeModel::GetUnpairedElectronCoordinates( *sl));
        BCL_MessageDbg( "calculated unpaired electron coordinates " + util::Format()( result));
        const linal::Vector3D expected( 39.914, -2.778, -0.139);
        BCL_MessageDbg( "expected unpaired electron coordinates " + util::Format()( expected));
        BCL_ExampleCheckWithinTolerance( result, expected, 0.0001);
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintConeModel

  const ExampleClass::EnumType ExampleRestraintConeModel::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintConeModel())
  );

} // namespace bcl
