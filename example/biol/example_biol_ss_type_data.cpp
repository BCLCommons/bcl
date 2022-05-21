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
#include "biol/bcl_biol_ss_type_data.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_ss_type_data.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolSSTypeData :
    public ExampleInterface
  {
  public:

    ExampleBiolSSTypeData *Clone() const
    {
      return new ExampleBiolSSTypeData( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::SSTypeData sstypedata_default;

      const char one_letter_code( 'H');
      const bool is_structured( true);
      const double radial_extent( 4.240);
      const double angle_per_turn( -100.0 / 180.0 * math::g_Pi);
      const double rise_z_per_residue( 1.50247);
      const double phi( -0.996767);
      const double psi( -0.80718);
      const size_t fragment( 5);
      const size_t window( 4);
      const linal::Vector3D three_state_prediction( 1.0, 0.0, 0.0);
      const math::Range< double> phi_range( math::Angle::Radian( -135.0), math::Angle::Radian( -25.0));
      const math::Range< double> psi_range( math::Angle::Radian( -70.0), math::Angle::Radian( 20.0));

      // construct from information
      biol::SSTypeData sstypedata_constr
      (
        one_letter_code,
        is_structured,
        radial_extent,
        angle_per_turn,
        rise_z_per_residue,
        phi,
        psi,
        fragment,
        window,
        three_state_prediction,
        phi_range,
        psi_range
      );

      // copy
      biol::SSTypeData sstypedata_copy( sstypedata_constr);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( sstypedata_constr.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( ptr->GetClassIdentifier(), GetStaticClassName< biol::SSTypeData>());

      // one letter code
      BCL_ExampleCheck( sstypedata_constr.GetOneLetterCode(), one_letter_code);

      // is structured
      BCL_MessageStd( "sstype is structured: " + util::Format()( sstypedata_constr.IsStructured()));
      BCL_Example_Check
      (
        sstypedata_constr.IsStructured() == is_structured,
        "is structured incorrectly initialized"
      );

      // radial extent
      BCL_MessageStd( "sstype radial extent: " + util::Format()( sstypedata_constr.GetRadialExtent()));
      BCL_Example_Check
      (
        sstypedata_constr.GetRadialExtent() == radial_extent,
        "radial extent incorrectly initialized"
      );

      // angle per turn
      BCL_MessageStd( "sstype angle per turn: " + util::Format()( sstypedata_constr.GetAnglePerTurn()));
      BCL_Example_Check
      (
        sstypedata_constr.GetAnglePerTurn() == angle_per_turn,
        "angle per turn incorrectly initialized"
      );

      // angle per turn
      BCL_MessageStd( "sstype rise in z per residue: " + util::Format()( sstypedata_constr.GetRiseInZPerResidue()));
      BCL_Example_Check
      (
        sstypedata_constr.GetRiseInZPerResidue() == rise_z_per_residue,
        "rise in z incorrectly initialized"
      );

      // ideal phi angle
      BCL_MessageStd( "sstype ideal phi: " + util::Format()( sstypedata_constr.GetIdealPhi()));
      BCL_Example_Check
      (
        sstypedata_constr.GetIdealPhi() == phi,
        "ideal phi incorrectly initialized"
      );

      // ideal psi angle
      BCL_MessageStd( "sstype ideal psi: " + util::Format()( sstypedata_constr.GetIdealPsi()));
      BCL_Example_Check
      (
        sstypedata_constr.GetIdealPsi() == psi,
        "ideal psi incorrectly initialized"
      );

      // fragment length
      BCL_MessageStd( "sstype fragment length: " + util::Format()( sstypedata_constr.GetFragmentLength()));
      BCL_Example_Check
      (
        sstypedata_constr.GetFragmentLength() == fragment,
        "fragment length incorrectly initialized"
      );

      // contact window radius
      BCL_MessageStd( "sstype contact window radius: " + util::Format()( sstypedata_constr.GetContactWindowRadius()));
      BCL_Example_Check
      (
        sstypedata_constr.GetContactWindowRadius() == window,
        "contact window radius incorrectly initialized"
      );

      // contact window length
      BCL_MessageStd( "sstype contact window length: " + util::Format()( sstypedata_constr.GetContactWindowLength()));
      BCL_Example_Check
      (
        sstypedata_constr.GetContactWindowLength() == 2 * window + 1,
        "contact window length incorrect"
      );

      const math::TransformationMatrix3D transformation
      (
        math::TransformationMatrix3D( linal::Vector3D( 0.0, 0.0, rise_z_per_residue))
        ( math::RotationMatrix3D( coord::GetAxes().e_Z, angle_per_turn))
      );
      // transformation matrix for residue
      BCL_MessageStd( "sstype residue transformation matrix: " + util::Format()( sstypedata_constr.GetTransformationMatrixForResidues()));
      BCL_Example_Check
      (
        sstypedata_constr.GetTransformationMatrixForResidues() == transformation,
        "transformation incorrect"
      );

      // check three state vector
      BCL_ExampleCheck( sstypedata_constr.GetThreeStatePrediction(), three_state_prediction);

      // check range
      BCL_ExampleCheck( sstypedata_constr.GetBackbonePhiRange(), phi_range);
      BCL_ExampleCheck( sstypedata_constr.GetBackbonePsiRange(), psi_range);

      // check default
      // compare written and read object
      BCL_Example_Check
      (
        !sstypedata_default.IsStructured()                             &&
        !util::IsDefined( sstypedata_default.GetAnglePerTurn())        &&
        !util::IsDefined( sstypedata_default.GetContactWindowRadius()) &&
        !util::IsDefined( sstypedata_default.GetFragmentLength())      &&
        !util::IsDefined( sstypedata_default.GetIdealPhi())            &&
        !util::IsDefined( sstypedata_default.GetIdealPsi())            &&
        sstypedata_default.GetOneLetterCode() == ' '                   &&
        !util::IsDefined( sstypedata_default.GetRadialExtent())        &&
        !util::IsDefined( sstypedata_default.GetRiseInZPerResidue())   &&
        sstypedata_default.GetTransformationMatrixForResidues() == math::TransformationMatrix3D() &&
        !sstypedata_default.GetThreeStatePrediction().IsDefined() &&
        sstypedata_default.GetBackbonePhiRange().IsEmpty() &&
        sstypedata_default.GetBackbonePsiRange().IsEmpty(),
        "default sstypedata is incorrect"
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

      BCL_MessageStd( "testing read and write functionalities for biol::SSTypeData");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( sstypedata_constr);
      BCL_MessageVrb( "read object");
      biol::SSTypeData sstypedata_read;
      ReadBCLObject( sstypedata_read);

      // compare written and read object
      BCL_Example_Check
      (
        sstypedata_constr.IsStructured()                       == sstypedata_read.IsStructured()            &&
        math::EqualWithinTolerance( sstypedata_constr.GetAnglePerTurn(), sstypedata_read.GetAnglePerTurn()) &&
        sstypedata_constr.GetContactWindowLength()             == sstypedata_read.GetContactWindowLength()  &&
        sstypedata_constr.GetContactWindowRadius()             == sstypedata_read.GetContactWindowRadius()  &&
        sstypedata_constr.GetFragmentLength()                  == sstypedata_read.GetFragmentLength()       &&
        sstypedata_constr.GetIdealPhi()                        == sstypedata_read.GetIdealPhi()             &&
        sstypedata_constr.GetIdealPsi()                        == sstypedata_read.GetIdealPsi()             &&
        sstypedata_constr.GetOneLetterCode()                   == sstypedata_read.GetOneLetterCode()        &&
        sstypedata_constr.GetRadialExtent()                    == sstypedata_read.GetRadialExtent()         &&
        sstypedata_constr.GetRiseInZPerResidue()               == sstypedata_read.GetRiseInZPerResidue()    &&
        math::EqualWithinTolerance
        (
          linal::Vector< double>( 16, sstypedata_constr.GetTransformationMatrixForResidues().GetMatrix().Begin()),
          linal::Vector< double>( 16, sstypedata_read.GetTransformationMatrixForResidues().GetMatrix().Begin())
        ) &&
        sstypedata_constr.GetThreeStatePrediction()            == sstypedata_read.GetThreeStatePrediction() &&
        sstypedata_constr.GetBackbonePhiRange().GetString()    == sstypedata_read.GetBackbonePhiRange().GetString() &&
        sstypedata_constr.GetBackbonePsiRange().GetString()    == sstypedata_read.GetBackbonePsiRange().GetString(),
        "read sstypedata is different from written: " + util::Format()( sstypedata_constr) + " != \n" +
        util::Format()( sstypedata_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolSSTypeData

  const ExampleClass::EnumType ExampleBiolSSTypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolSSTypeData())
  );

} // namespace bcl

