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
#include "sdf/bcl_sdf_atom_info.h"

// includes from bcl - sorted alphabetically
#include "sdf/bcl_sdf_mdl_line_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_atom_info.cpp
  //!
  //! @author mendenjl
  //! @date Mar 02, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfAtomInfo :
    public ExampleInterface
  {
  public:

    ExampleSdfAtomInfo *Clone() const
    {
      return new ExampleSdfAtomInfo( *this);
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
      sdf::AtomInfo carbon_te_s( chemistry::GetAtomTypes().C_TeTeTeTe, chemistry::e_SChirality);
      linal::Vector3D test_coordinates( 1.2, 2.4, 4.8);
      sdf::AtomInfo carbon_tr_with_position_no_h
      (
        chemistry::GetAtomTypes().C_TrTrPi,
        chemistry::e_NonChiral,
        test_coordinates,
        false
      );

    /////////////////
    // data access //
    /////////////////

      // check default constructor
      BCL_ExampleCheck( sdf::AtomInfo().GetAtomType(),    chemistry::AtomType());
      BCL_ExampleCheck( sdf::AtomInfo().GetChirality(),   chemistry::e_UnknownChirality);
      BCL_ExampleCheck( sdf::AtomInfo().GetCoordinates(), linal::Vector3D( 0.0));
      BCL_ExampleCheck( sdf::AtomInfo().CanAddH(),        true);
      //BCL_ExampleCheck( sdf::AtomInfo().GetAtomMapping(), 0);
      BCL_ExampleCheck( sdf::AtomInfo().ToMdlAtomLine(),  sdf::GetDefaultLine( sdf::e_AtomLine));

      // check constructor from atom type, chirality
      BCL_ExampleCheck( carbon_te_s.GetAtomType(),    chemistry::GetAtomTypes().C_TeTeTeTe);
      BCL_ExampleCheck( carbon_te_s.GetChirality(),   chemistry::e_SChirality);
      BCL_ExampleCheck( carbon_te_s.CanAddH(),        true);
      BCL_ExampleCheck( carbon_tr_with_position_no_h.GetCoordinates(), test_coordinates);
      BCL_ExampleCheck( carbon_tr_with_position_no_h.CanAddH(),        false);

      // test get string on a non-trivial example
      const std::string expected_carbon_te_s("    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0");
      const std::string expected_carbon_tr_with_position
      (
        "    1.2000    2.4000    4.8000 C   0  3  0  1  0  0  0  0  0  0  0  0"
      );
      BCL_ExampleCheck( carbon_te_s.ToMdlAtomLine(), expected_carbon_te_s);
      BCL_ExampleCheck( carbon_tr_with_position_no_h.ToMdlAtomLine(), expected_carbon_tr_with_position);
      // test constructor from string; output should differ in valence, since with only the element type and
      // no mdl bonds, there is no way to determine the valence
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_te_s).ToMdlAtomLine(),
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
      );
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_te_s).GetAtomType()->GetElementType(),
        chemistry::GetElementTypes().e_Carbon
      );
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position).ToMdlAtomLine(),
        "    1.2000    2.4000    4.8000 C   0  3  0  1  0  0  0  0  0  0  0  0"
      );
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position).GetAtomType()->GetElementType(),
        chemistry::GetElementTypes().e_Carbon
      );
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position).GetAtomType()->GetFormalCharge(),
        1
      );
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position).GetCoordinates(),
        test_coordinates
      );
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position).CanAddH(),
        false
      );
      /*
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position).GetAtomMapping(),
        0
      );
      */

      // test non-zero atom mapping
      sdf::AtomInfo carbon_tr_with_position_no_h_map
      (
        chemistry::GetAtomTypes().C_TrTrPi,
        chemistry::e_NonChiral,
        test_coordinates,
        false
        //1 // mapping
      );

      //BCL_ExampleCheck( carbon_tr_with_position_no_h_map.GetAtomMapping(), 1);
      //carbon_tr_with_position_no_h_map.SetAtomMapping( 2);

      /*
      const std::string expected_carbon_tr_with_position_map_line
      (
        "    1.2000    2.4000    4.8000 C   0  3  0  1  0  0  0  0  0  2  0  0"
      );
      */

      /*
      BCL_ExampleCheck
      (
        sdf::AtomInfo().ExtractMdlAtomLineInfo( expected_carbon_tr_with_position_map_line).GetAtomMapping(),
        carbon_tr_with_position_no_h_map.GetAtomMapping()
      );
      */
      /*
      BCL_ExampleCheck
      (
        carbon_tr_with_position_no_h_map.ToMdlAtomLine(),
        expected_carbon_tr_with_position_map_line
      );
      */

    //////////////////////
    // input and output //
    //////////////////////

      // test io
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( carbon_tr_with_position_no_h, sdf::AtomInfo()), true);
      BCL_ExampleCheck( TestBCLObjectOutputDiffers( carbon_tr_with_position_no_h, sdf::AtomInfo()), true);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfAtomInfo

  const ExampleClass::EnumType ExampleSdfAtomInfo::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfAtomInfo())
  );

} // namespace bcl
