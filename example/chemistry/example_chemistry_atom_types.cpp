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
#include "chemistry/bcl_chemistry_atom_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_types.cpp
  //!
  //! @author mendenjl
  //! @date July 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomTypes :
    public ExampleInterface
  {
  public:

    ExampleChemistryAtomTypes *Clone() const
    {
      return new ExampleChemistryAtomTypes( *this);
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
      chemistry::GetAtomTypes();

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier
      BCL_ExampleCheck( chemistry::GetAtomTypes().GetClassIdentifier(), "bcl::chemistry::AtomTypes");

      // test ability to dynamically add to the atom types by calling get atom type with something that is not a
      // defined gasteiger type
      BCL_ExampleCheck
      (
        chemistry::GetAtomTypes().GetAtomType( chemistry::GetElementTypes().e_Calcium, 4)->GetFormalCharge(),
        4
      );
      BCL_ExampleCheck
      (
        chemistry::GetAtomTypes().GetAtomType( chemistry::GetElementTypes().e_Calcium, 4)->GetElementType(),
        chemistry::GetElementTypes().e_Calcium
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

      // write the output from the AtomTypeData public functions for all atom types
      for
      (
        chemistry::AtomTypes::const_iterator
          itr( chemistry::GetAtomTypes().Begin()),
          itr_end( chemistry::GetAtomTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        chemistry::AtomType type( *itr);
        const chemistry::AtomTypeData &data( **itr);
        BCL_MessageDbg
        (
          "AtomType " + type.GetName()
          + " " + data.GetElementType().GetName()
          + " " + util::Format()( data.GetHybridOrbitalType().GetIndex())
          + " " + util::Format()( data.GetFormalCharge())
          + " " + util::Format()( data.GetNumberBonds())
          + " " + util::Format()( data.GetNumberHybridBonds())
          + " " + util::Format()( data.GetNumberElectronsInBonds())
          + " " + util::Format()( data.GetNumberHybridOrbitals())
          + " " + util::Format()( data.GetNumberHybridLonePairs())
          + " " + util::Format()( data.GetNumberElectronsInPOrbitals())
          + " " + util::Format()( data.GetNumberPiOrbitals())
          + " " + util::Format()( data.GetNumberSigmaOrbitals())
          + " " + util::Format()( data.GetNumberUnhybridizedSigmaOrbitals())
          + " " + util::Format()( data.GetNumberUnhybridizedLonePairs())
          + " " + util::Format()( data.GetValenceElectronsSP())
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleChemistryAtomTypes

  const ExampleClass::EnumType ExampleChemistryAtomTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomTypes())
  );
} // namespace bcl
