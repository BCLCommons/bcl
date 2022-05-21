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
#include "nmr/bcl_nmr_signal_1d.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_signal_1d.cpp
  //!
  //! @author mueller
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrSignal1D :
    public ExampleInterface
  {
  public:

    ExampleNmrSignal1D *Clone() const
    {
      return new ExampleNmrSignal1D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      storage::Vector< sdf::AtomInfo> complete_intializer
      (
        storage::Vector< sdf::AtomInfo>::Create
        (
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_NonChiral, linal::Vector3D( double( 0.0))),
          sdf::AtomInfo( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_NonChiral, linal::Vector3D( double( -1.0)))
        )
      );

      // construct bond connectivities
      storage::Vector< sdf::BondInfo> atom_bonds
      (
          1, sdf::BondInfo( 0, 1, chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond)
      );

      chemistry::AtomVector< chemistry::AtomComplete> atom_vector( complete_intializer, atom_bonds);
      storage::Pair< nmr::Signal1D::PatternTypeEnum, double> pattern_type_coupling( nmr::Signal1D::e_Doublet, 5);
      const storage::Vector< storage::Pair< nmr::Signal1D::PatternTypeEnum, double> > pattern( 1, pattern_type_coupling);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test standard constructor
      nmr::Signal1D signal_1d_a;

      nmr::Signal1D signal_1d_b
      (
        120.0,
        1.0,
        pattern,
        atom_vector( 0)
      );

      BCL_MessageDbg
      (
        "Signal1D: "
        + util::Format()( signal_1d_b)
      );

    /////////////////
    // data access //
    /////////////////

      signal_1d_a.SetChemicalShift( 100);

      BCL_Example_Check
      (
        signal_1d_a.GetChemicalShift() == 100 && signal_1d_b.GetChemicalShift() == 120,
        "GetChemicalShift() gives wrong result."
      );

      signal_1d_a.SetIntegral( 3);

      BCL_Example_Check
      (
        signal_1d_a.GetIntegral() == 3 && signal_1d_b.GetIntegral() == 1,
        "GetIntegral() gives wrong result."
      );

      signal_1d_a.SetAtomInvolvedInSignal( atom_vector( 0));

      BCL_Example_Check
      (
        signal_1d_a.GetAtomInvolvedInSignal() == util::ToSiPtr( atom_vector( 0))
        && signal_1d_b.GetAtomInvolvedInSignal() == util::ToSiPtr( atom_vector( 0)),
        "GetAtomInvolvedInSignal() gives wrong result."
      );

      signal_1d_a.SetPattern( pattern);

      BCL_Example_Check
      (
        signal_1d_a.GetPattern()( 0).First() == nmr::Signal1D::e_Doublet
        && signal_1d_a.GetPattern()( 0).Second() == 5,
        "GetPattern() gives wrong result."
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

      // write object
      WriteBCLObject( signal_1d_a);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrSignal1D

  const ExampleClass::EnumType ExampleNmrSignal1D::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrSignal1D())
  );

} // namespace bcl
