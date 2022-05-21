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
#include "quality/bcl_quality_dme.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_dme.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityDME :
    public ExampleInterface
  {
  public:

    ExampleQualityDME *Clone() const
    {
      return new ExampleQualityDME( *this);
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
      storage::Set< biol::AtomType> atom_types( biol::GetAtomTypes().CA);
      util::SiPtrVector< const linal::Vector3D> native_coords( native.GetAtomCoordinates( atom_types));
      util::SiPtrVector< const linal::Vector3D> model_coords( model.GetAtomCoordinates( atom_types));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      quality::DME dme;

      // clone
      util::ShPtr< quality::DME> sp_dme( dme.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_ExampleCheck( GetStaticClassName< quality::DME>(), dme.GetClassIdentifier());

    ////////////////
    // operations //
    ////////////////

      BCL_MessageStd( "test operator()");
      // expected value
      const double expected_dme( 2.31551);
      const double dme_value( dme.CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "dme value: " + util::Format()( dme_value));
      BCL_ExampleCheckWithinTolerance( dme_value, expected_dme, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityDME

  const ExampleClass::EnumType ExampleQualityDME::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityDME())
  );
  
} // namespace bcl
