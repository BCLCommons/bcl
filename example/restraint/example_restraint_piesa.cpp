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
#include "restraint/bcl_restraint_piesa.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_piesa.cpp
  //! @brief this example tests the implementation of the class scoring simulating PIESA spectra from protein models
  //!
  //! @author fischea
  //! @date Nov 22, 2015
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintPiesa :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to a new ExampleRestraintPiesa
    //! @return pointer to a new ExampleRestraintPiesa
    ExampleRestraintPiesa *Clone() const
    {
      return new ExampleRestraintPiesa( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a PIESA object
      restraint::Piesa piesa;

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( piesa.GetClassIdentifier(), GetStaticClassName< restraint::Piesa>());

    ////////////////
    // operations //
    ////////////////

      // create a protein model for testing
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Restraint, "1py6_idealized.pdb"));
      assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDBFilename( pdb_filename));
      biol::AABackBoneCompleter aa_bb_compl( true, false, false);
      model = *aa_bb_compl.CompleteProteinModel( model);
      util::SiPtrVector< const biol::AABase> aas( model.GetAminoAcids());

      // simulate the PIESA spectrum for SSEs in the protein model
      util::SiPtrVector< const assemble::SSE> sses( model.GetSSEs( biol::GetSSTypes().HELIX));
      util::ShPtr< assemble::SSE> sp_sse( util::CloneToShPtr( *sses( 0)));
      sp_sse->SetOrigin( linal::Vector3D( 0.0, 0.0, 0.0));

      linal::Vector3D axis( sp_sse->GetMainAxis().GetDirection());
      // axis = linal::Vector3D( 0.0, -10.0, 0.0);
      axis.Normalize();

      linal::Vector3D tmp( 1.0, 0.0, 0.0);
      tmp = ( tmp - ( tmp * axis) * axis).Normalize();

      sp_sse->Rotate( sp_sse->GetCenter(), tmp, 90.0 / 180.0 * math::g_Pi);

      util::ShPtr< restraint::Piesa> sp_piesa( restraint::Piesa::Create( *sp_sse, axis));

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; // class ExampleRestraintPiesa

  //! single instance of this class
  const ExampleClass::EnumType ExampleRestraintPiesa::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintPiesa())
  );

} // namespace bcl
