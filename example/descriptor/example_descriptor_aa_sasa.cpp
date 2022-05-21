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
#include "descriptor/bcl_descriptor_aa_sasa.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_aa_sasa.cpp
  //!
  //! @author woetzen, mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorAASasa :
    public ExampleInterface
  {
  public:

    ExampleDescriptorAASasa *Clone() const
    {
      return new ExampleDescriptorAASasa( *this);
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

      // default constructor
      descriptor::AASasa aa_sasa( assemble::AASasaOLS(), "", "");
      aa_sasa.SetObject( protein_model);

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( aa_sasa.GetClassIdentifier(), GetStaticClassName( aa_sasa));

      // test the GetLength function
      BCL_ExampleCheck( aa_sasa.GetSizeOfFeatures(), 1);
      BCL_ExampleCheck( aa_sasa.GetType().GetDimension(), 1);

    ///////////////
    // operators //
    ///////////////

      // test the operator
      BCL_MessageStd
      (
        "exposure for amino acid: " + itr_amino_acid( 0)->GetIdentification() + " = "
        + util::Format()( aa_sasa( itr_amino_acid))
      );

      const double expected_aa_sasa( 0.06);
      // check the descriptions produced
      BCL_ExampleCheckWithinTolerance( aa_sasa( itr_amino_acid)( 0), expected_aa_sasa, 0.0001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorAASasa

  const ExampleClass::EnumType ExampleDescriptorAASasa::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorAASasa())
  );

} // namespace bcl
