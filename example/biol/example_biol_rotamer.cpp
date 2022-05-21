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
#include "biol/bcl_biol_rotamer.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_rotamer.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Aug 25, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolRotamer :
    public ExampleInterface
  {
  public:

    ExampleBiolRotamer *Clone() const
    {
      return new ExampleBiolRotamer( *this);
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
      biol::Rotamer default_constr;
      BCL_ExampleCheck( default_constr.IsEmpty(), true);

      // clone constructor
      default_constr.Insert( biol::ChiAngle( biol::ChiAngle::e_Four, 45.0, math::Angle::e_Degree));
      util::ShPtr< biol::Rotamer> clone_constr( default_constr.Clone());
      BCL_ExampleCheck( default_constr.IsEmpty(), false);

    /////////////////
    // data access //
    /////////////////

      // IsEmpty
      BCL_ExampleCheck( default_constr.IsEmpty(), false);

      // Begin
      BCL_ExampleCheck( default_constr.Begin()->GetChi(), biol::ChiAngle::e_Four);

      // End
      BCL_ExampleCheck( ( --default_constr.End())->GetChi(), biol::ChiAngle::e_Four);

      // GetSize
      BCL_ExampleCheck( default_constr.GetSize(), 1);

      // GetChis
      BCL_ExampleCheck( default_constr.GetChis().GetSize(), 1);
      BCL_ExampleCheck( *default_constr.GetChis().Begin(), biol::ChiAngle::e_Four);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // Insert
      {
        std::pair< biol::Rotamer::const_iterator, bool> success
        (
          default_constr.Insert( biol::ChiAngle( biol::ChiAngle::e_Four, 45.0, math::Angle::e_Degree))
        );

        BCL_ExampleCheck( success.second, false);
      }
      {
        std::pair< biol::Rotamer::const_iterator, bool> success
        (
          default_constr.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, 90.0, math::Angle::e_Degree))
        );

        BCL_ExampleCheck( success.second, true);
        BCL_ExampleCheck( success.first->GetChi(), biol::ChiAngle::e_Two);
      }

      // GetAngle
      BCL_ExampleCheck( default_constr.GetAngle( biol::ChiAngle::e_Two, math::Angle::e_Degree), 90.0);
      BCL_ExampleCheck
      (
        util::IsDefined( default_constr.GetAngle( biol::ChiAngle::e_One, math::Angle::e_Degree)), false
      );

      // ChiMatchDependent
      {
        biol::Rotamer rotamer_a;
        rotamer_a.Insert( biol::ChiAngle( biol::ChiAngle::e_One, 90.0, math::Angle::e_Degree));
        rotamer_a.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, 80.0, math::Angle::e_Degree));
        biol::Rotamer rotamer_b;
        rotamer_b.Insert( biol::ChiAngle( biol::ChiAngle::e_One, 90.0, math::Angle::e_Degree));
        rotamer_b.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, -120.0, math::Angle::e_Degree));
        const storage::Set< biol::ChiAngle::ChiEnum> chi
        (
          rotamer_a.ChiMatchDependent( rotamer_b, math::Angle::e_Degree, 5)
        );
        BCL_ExampleCheck( chi.GetSize(), 1);
        BCL_ExampleCheck( *chi.Begin(), biol::ChiAngle::e_One);
      }
      {
        biol::Rotamer rotamer_a;
        rotamer_a.Insert( biol::ChiAngle( biol::ChiAngle::e_One, 90.0, math::Angle::e_Degree));
        rotamer_a.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, 80.0, math::Angle::e_Degree));
        biol::Rotamer rotamer_b;
        rotamer_b.Insert( biol::ChiAngle( biol::ChiAngle::e_One, 130.0, math::Angle::e_Degree));
        rotamer_b.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, 81.0, math::Angle::e_Degree));
        const storage::Set< biol::ChiAngle::ChiEnum> chi
        (
          rotamer_a.ChiMatchDependent( rotamer_b, math::Angle::e_Degree, 5)
        );
        BCL_ExampleCheck( chi.IsEmpty(), true);
      }
      {
        biol::Rotamer rotamer_a;
        rotamer_a.Insert( biol::ChiAngle( biol::ChiAngle::e_One, 90.0, math::Angle::e_Degree));
        rotamer_a.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, 80.0, math::Angle::e_Degree));
        biol::Rotamer rotamer_b;
        rotamer_b.Insert( biol::ChiAngle( biol::ChiAngle::e_One, 92.0, math::Angle::e_Degree));
        rotamer_b.Insert( biol::ChiAngle( biol::ChiAngle::e_Two, 81.0, math::Angle::e_Degree));
        const storage::Set< biol::ChiAngle::ChiEnum> chi
        (
          rotamer_a.ChiMatchDependent( rotamer_b, math::Angle::e_Degree, 5)
        );
        BCL_ExampleCheck( chi.GetSize(), 2);
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( default_constr);
      biol::Rotamer read_rotamer;
      ReadBCLObject( read_rotamer);
      BCL_ExampleCheck( read_rotamer.IsEmpty(), false);
      BCL_ExampleCheck( default_constr.GetAngle( biol::ChiAngle::e_Two, math::Angle::e_Degree), 90.0);
      BCL_ExampleCheck( default_constr.GetAngle( biol::ChiAngle::e_Four, math::Angle::e_Degree), 45.0);

    //////////////////////
    // helper functions //
    //////////////////////

      // ReadSimple
      {
        const std::string input_filename
        (
          AddExampleInputPathToFilename( e_Biology, "bcl_biol_rotamer_read_simple.txt")
        );
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, input_filename);
        biol::Rotamer rotamer;
        rotamer.ReadSimple( read);
        BCL_ExampleCheck( rotamer.GetSize(), 2);
        BCL_ExampleCheck( rotamer.GetAngle( biol::ChiAngle::e_Two, math::Angle::e_Degree), 90);
        BCL_ExampleCheck( rotamer.GetAngle( biol::ChiAngle::e_Four, math::Angle::e_Degree), 45);
      }

      // ChiAngleLessThan
      BCL_ExampleCheck
      (
        biol::Rotamer::ChiAngleLessThan()
        (
          biol::ChiAngle( biol::ChiAngle::e_Four, 45.0, math::Angle::e_Degree),
          biol::ChiAngle( biol::ChiAngle::e_Four, 45.0, math::Angle::e_Degree)
        ),
        false
      );
      BCL_ExampleCheck
      (
        biol::Rotamer::ChiAngleLessThan()
        (
          biol::ChiAngle( biol::ChiAngle::e_Five, 45.0, math::Angle::e_Degree),
          biol::ChiAngle( biol::ChiAngle::e_Four, 45.0, math::Angle::e_Degree)
        ),
        false
      );
      BCL_ExampleCheck
      (
        biol::Rotamer::ChiAngleLessThan()
        (
          biol::ChiAngle( biol::ChiAngle::e_Three, 45.0, math::Angle::e_Degree),
          biol::ChiAngle( biol::ChiAngle::e_Four, 45.0, math::Angle::e_Degree)
        ),
        true
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolRotamer

  const ExampleClass::EnumType ExampleBiolRotamer::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolRotamer())
  );

} // namespace bcl
