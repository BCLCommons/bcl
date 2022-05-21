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
#include "descriptor/bcl_descriptor_aa_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_phi_psi.cpp
  //!
  //! @author mendenjl
  //! @date Apr 23, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAAPhiPsi :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAAPhiPsi *Clone() const
    {
      return new ExampleDescriptorAAPhiPsi( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      assemble::ProteinModelWithCache protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone),
        true
      );

      // create an iterator on the protein model
      descriptor::Iterator< biol::AABase> itr_amino_acid
      (
        descriptor::Type( 1, false, descriptor::Type::e_Symmetric),
        protein_model
      );

      // goto the 21st position in the protein
      itr_amino_acid.GotoPosition( 20);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor for phi and psi
      descriptor::AAPhiPsi aa_phi( true, -180.0), aa_psi( false, -90.0);
      aa_phi.SetObject( protein_model);
      aa_psi.SetObject( protein_model);

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( aa_phi.GetClassIdentifier(), GetStaticClassName( aa_phi));

      // test the basic interface functions
      BCL_ExampleCheck( aa_psi.GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( aa_phi.GetType().GetDimension(), 1);

    ///////////////
    // operators //
    ///////////////

      const double expected_aa_phi( -74.211), expected_aa_psi( 143.658);
      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance( aa_phi( itr_amino_acid)( 0), expected_aa_phi, 0.01);
      BCL_ExampleCheckWithinTolerance( aa_psi( itr_amino_acid)( 0), expected_aa_psi, 0.01);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAAPhiPsi

  const ExampleClass::EnumType ExampleDescriptorAAPhiPsi::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAAPhiPsi())
  );

} // namespace bcl
