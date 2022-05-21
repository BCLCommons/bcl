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
#include "chemistry/bcl_chemistry_atom_type_data.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_type_data.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomTypeData :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomTypeData *Clone() const
    {
      return new ExampleChemistryAtomTypeData( *this);
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

      // test default constructor
      chemistry::AtomTypeData atom_type_data_default;

      // test constructor given geometry and hybridization
      chemistry::AtomTypeData atom_type_data_chlorine
      (
        chemistry::GetElementTypes().e_Chlorine,
        chemistry::GetHybridOrbitalTypes().e_Unhybridized,
        0,
        0,
        storage::Set< chemistry::AtomicOrbitalTypesEnum>( chemistry::e_Pz),
        storage::Set< chemistry::AtomicOrbitalTypesEnum>::Create
        (
          chemistry::e_S,
          chemistry::e_Px,
          chemistry::e_Py
        ),
        15.03,
        3.73,
        util::GetUndefined< double>(),
        util::GetUndefined< double>(),
        util::GetUndefined< double>(),
        util::GetUndefined< double>(),
        util::GetUndefined< double>()
      );

      // test copy constructor
      chemistry::AtomTypeData atom_type_data_chlorine_copy( atom_type_data_chlorine);

      // test Clone
      util::ShPtr< chemistry::AtomTypeData> sp_atom_type_data_chlorine_clone( atom_type_data_chlorine.Clone());

    /////////////////
    // data access //
    /////////////////

      // test GetElementType
      BCL_ExampleCheck( atom_type_data_chlorine.GetElementType(), chemistry::GetElementTypes().e_Chlorine);

      // test GetNumberHybridBonds
      BCL_ExampleCheck( atom_type_data_chlorine.GetNumberHybridBonds(), 1);

      // test GetNumberBonds
      BCL_ExampleCheck( atom_type_data_chlorine.GetNumberBonds(), 1);

      // test GetNumberPiOrbitals
      BCL_ExampleCheck( atom_type_data_chlorine.GetNumberPiOrbitals(), 1);

      // test GetFormalCharge
      BCL_ExampleCheck( atom_type_data_chlorine.GetFormalCharge(), 0);

      // test GetNumberElectronsInBonds
      BCL_ExampleCheck( atom_type_data_chlorine.GetNumberElectronsInBonds(), 1);

      // test that the electronegativity values were set correctly
      BCL_ExampleCheckWithinAbsTolerance
      (
        atom_type_data_chlorine.GetAtomTypeProperty( chemistry::AtomTypeData::e_SigmaOrbitalElectronegativityPauling),
        2.94504,
        0.00001
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test writing atom_type_data to file
      WriteBCLObject( atom_type_data_chlorine);
      // test reading atom_type_data from file
      chemistry::AtomTypeData read_atom_type_data_chlorine;
      ReadBCLObject( read_atom_type_data_chlorine);

      BCL_ExampleCheck( read_atom_type_data_chlorine.GetElementType(), chemistry::GetElementTypes().e_Chlorine);
      BCL_ExampleCheck( read_atom_type_data_chlorine.GetFormalCharge(), 0);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryAtomTypeData

  const ExampleClass::EnumType ExampleChemistryAtomTypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomTypeData())
  );

} // namespace bcl
