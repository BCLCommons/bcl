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
#include "sdf/bcl_sdf_bond_info.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_line_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_bond_info.cpp
  //!
  //! @author mendenjl
  //! @date Mar 02, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfBondInfo :
    public ExampleInterface
  {
  public:

    ExampleSdfBondInfo *Clone() const
    {
      return new ExampleSdfBondInfo( *this);
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

      // constructors from members
      sdf::BondInfo two_one_dbl_unk_iso( 2, 1,     chemistry::GetConstitutionalBondTypes().e_ConjugatedDoubleBond);
      sdf::BondInfo zero_four_dbl_e_iso( 0, 4,     chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E);

    /////////////////
    // data access //
    /////////////////

      // check default constructor
      BCL_ExampleCheck( sdf::BondInfo().ToMdlBondLine(), sdf::GetDefaultLine( sdf::e_BondLine));

      // check that constructors always order the atom indices correctly
      BCL_ExampleCheck( two_one_dbl_unk_iso.GetAtomIndexHigh(), 2);
      BCL_ExampleCheck( two_one_dbl_unk_iso.GetAtomIndexLow(),  1);
      BCL_ExampleCheck( zero_four_dbl_e_iso.GetAtomIndexHigh(), 4);
      BCL_ExampleCheck( zero_four_dbl_e_iso.GetAtomIndexLow(),  0);

      // test that bond types are accurate
      BCL_ExampleCheck
      (
        zero_four_dbl_e_iso.GetConfigurationalBondType(),
        chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E
      );
      BCL_ExampleCheck
      (
        two_one_dbl_unk_iso.GetConstitutionalBondType(),
        chemistry::GetConstitutionalBondTypes().e_ConjugatedDoubleBond
      );

      // test is valid
      BCL_ExampleCheck( sdf::BondInfo().IsValid(), false);
      BCL_ExampleCheck( zero_four_dbl_e_iso.IsValid(), true);

    /////////////////
    // data access //
    /////////////////

      sdf::BondInfo two_one_dbl_non_iso( 2, 1,     chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond);
      sdf::BondInfo two_six_single(      2, 6,     chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBond);
      sdf::BondInfo two_one_single_aro(  2, 1,     chemistry::GetConfigurationalBondTypes().e_AromaticSingleBond);
      sdf::BondInfo three_digit_indices( 123, 456, chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond);

      // test set isometry
      two_one_dbl_unk_iso.SetIsometry( chemistry::e_EIsometry);
      BCL_ExampleIndirectCheck
      (
        two_one_dbl_unk_iso.GetConfigurationalBondType(),
        chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E,
        "two_one_dbl_unk_iso.SetIsometry( chemistry::e_EIsometry)"
      );
      two_one_dbl_unk_iso.SetIsometry( chemistry::e_UnknownIsometry);
      BCL_ExampleIndirectCheck
      (
        two_one_dbl_unk_iso.GetConfigurationalBondType(),
        chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond_X,
        "two_one_dbl_unk_iso.SetIsometry( chemistry::e_UnknownIsometry)"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // writing it out using GetString() should only give us the atom indices and bond orders;
      // isometry is stored elsewhere;  atom indices are offset by 1
      const std::string expected_two_one_dbl_non_iso_to_mdl_line( "  2  3  2  0  0  0  0");
      const std::string expected_two_six_single_to_mdl_line(      "  3  7  1  0  0  0  0");
      const std::string expected_two_one_single_aro_to_mdl_line(  "  2  3  1  0  0  0  0");
      const std::string expected_three_digit_indices_to_mdl_line( "124457  1  0  0  0  0");

      BCL_ExampleCheck( two_one_dbl_non_iso.ToMdlBondLine(), expected_two_one_dbl_non_iso_to_mdl_line);
      BCL_ExampleCheck( two_six_single.ToMdlBondLine(),      expected_two_six_single_to_mdl_line);
      BCL_ExampleCheck( two_one_single_aro.ToMdlBondLine(),  expected_two_one_single_aro_to_mdl_line);
      BCL_ExampleCheck( three_digit_indices.ToMdlBondLine(), expected_three_digit_indices_to_mdl_line);

      // test ExtractMdlBondLineInfo
      BCL_ExampleCheck
      (
        sdf::BondInfo().ExtractMdlBondLineInfo( expected_two_one_dbl_non_iso_to_mdl_line).ToMdlBondLine(),
        expected_two_one_dbl_non_iso_to_mdl_line
      );
      BCL_ExampleCheck
      (
        sdf::BondInfo().ExtractMdlBondLineInfo( expected_two_six_single_to_mdl_line).ToMdlBondLine(),
        expected_two_six_single_to_mdl_line
      );
      BCL_ExampleCheck
      (
        sdf::BondInfo().ExtractMdlBondLineInfo( expected_two_one_single_aro_to_mdl_line).ToMdlBondLine(),
        expected_two_one_single_aro_to_mdl_line
      );
      BCL_ExampleCheck
      (
        sdf::BondInfo().ExtractMdlBondLineInfo( expected_three_digit_indices_to_mdl_line).ToMdlBondLine(),
        expected_three_digit_indices_to_mdl_line
      );

      // test io
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( two_one_dbl_non_iso, sdf::BondInfo()), true);
      BCL_ExampleCheck( TestBCLObjectOutputDiffers( two_one_dbl_non_iso, sdf::BondInfo()), true);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfBondInfo

  const ExampleClass::EnumType ExampleSdfBondInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfBondInfo())
  );

} // namespace bcl
