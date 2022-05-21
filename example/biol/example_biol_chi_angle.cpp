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
#include "biol/bcl_biol_chi_angle.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_chi_angle.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Aug 25, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolChiAngle :
    public ExampleInterface
  {
  public:

    ExampleBiolChiAngle *Clone() const
    {
      return new ExampleBiolChiAngle( *this);
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

      // default constructor
      biol::ChiAngle def_constr;
      BCL_ExampleCheck( def_constr.GetChi(),   biol::ChiAngle::e_Undefined);
      BCL_ExampleCheck( util::IsDefined( def_constr.GetAngle( math::Angle::e_Radian)), false);

      // constructor taking parameter
      biol::ChiAngle param_constr( biol::ChiAngle::e_One);
      BCL_ExampleCheck( param_constr.GetChi(), biol::ChiAngle::e_One);
      BCL_ExampleCheck( util::IsDefined( param_constr.GetAngle( math::Angle::e_Radian)), false);

      // constructor taking parameters
      biol::ChiAngle param_constr_b( biol::ChiAngle::e_One, 180.0, math::Angle::e_Degree);
      BCL_ExampleCheck( param_constr_b.GetChi(), biol::ChiAngle::e_One);
      BCL_ExampleCheck( param_constr_b.GetAngle( math::Angle::e_Degree), 180.0);

      // clone constructor
      util::ShPtr< biol::ChiAngle> clone_constr( param_constr_b.Clone());
      BCL_ExampleCheck( clone_constr->GetChi(), biol::ChiAngle::e_One);
      BCL_ExampleCheck( clone_constr->GetAngle( math::Angle::e_Degree), 180.0);

    /////////////////
    // data access //
    /////////////////

      // GetChi
      BCL_ExampleCheck( param_constr_b.GetChi(), biol::ChiAngle::e_One);

      // GetAngle
      BCL_ExampleCheck( param_constr_b.GetAngle( math::Angle::e_Degree), 180.0);
      BCL_ExampleCheckWithinTolerance( param_constr_b.GetAngle( math::Angle::e_Radian), math::g_Pi, 0.001);
      {
        biol::ChiAngle test_get_angle( biol::ChiAngle::e_One, math::g_Pi / 2.0, math::Angle::e_Radian);
        BCL_ExampleCheck( test_get_angle.GetAngle( math::Angle::e_Degree), 90.0);
        BCL_ExampleCheck( test_get_angle.GetAngle( math::Angle::e_Radian), math::g_Pi / 2.0);
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // CalculateAngleDifference
      {
        biol::ChiAngle chi_a(  biol::ChiAngle::e_One, 45, math::Angle::e_Degree);
        biol::ChiAngle chi_b(  biol::ChiAngle::e_One, -30, math::Angle::e_Degree);
        BCL_ExampleCheck( chi_a.CalculateAngleDifference( chi_b, math::Angle::e_Degree), 75);
      }
      {
        biol::ChiAngle chi_a(  biol::ChiAngle::e_One, 45, math::Angle::e_Degree);
        biol::ChiAngle chi_b(  biol::ChiAngle::e_One, -30, math::Angle::e_Degree);
        const double calculated_difference( chi_a.CalculateAngleDifference( chi_b, math::Angle::e_Radian));
        const double radians( 1.3089969);
        BCL_MessageDbg
        (
          "calculated chi is " + util::Format()( calculated_difference) + " but expected " +
          util::Format()( radians)
        );
        BCL_ExampleCheckWithinTolerance( calculated_difference, radians, 0.001);
      }
      {
        biol::ChiAngle chi_a(  biol::ChiAngle::e_One, 45, math::Angle::e_Degree);
        biol::ChiAngle chi_b(  biol::ChiAngle::e_One, -160, math::Angle::e_Degree);
        BCL_ExampleCheck( chi_a.CalculateAngleDifference( chi_b, math::Angle::e_Degree), 155);
      }
      {
        biol::ChiAngle chi_a(  biol::ChiAngle::e_One, -175, math::Angle::e_Degree);
        biol::ChiAngle chi_b(  biol::ChiAngle::e_One, 160, math::Angle::e_Degree);
        BCL_ExampleCheck( chi_a.CalculateAngleDifference( chi_b, math::Angle::e_Degree), 25);
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( param_constr_b);
      biol::ChiAngle read_chi_angle;
      ReadBCLObject( read_chi_angle);
      BCL_ExampleCheck( read_chi_angle.GetChi(), biol::ChiAngle::e_One);
      BCL_ExampleCheck( read_chi_angle.GetAngle( math::Angle::e_Degree), 180.0);

    //////////////////////
    // helper functions //
    //////////////////////

      // ReadSimple
      {
        const std::string input_filename
        (
          AddExampleInputPathToFilename( e_Biology, "bcl_biol_chi_angle_read_simple.txt")
        );
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, input_filename);
        biol::ChiAngle chi;
        chi.ReadSimple( read);
        BCL_ExampleCheck( chi.GetChi(), biol::ChiAngle::e_Four);
        BCL_ExampleCheck( chi.GetAngle( math::Angle::e_Degree), 45);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolChiAngle

  const ExampleClass::EnumType ExampleBiolChiAngle::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolChiAngle())
  );

} // namespace bcl
