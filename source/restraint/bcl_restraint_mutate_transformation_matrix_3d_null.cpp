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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_null.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateTransformationMatrix3DNull::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateTransformationMatrix3DNull())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateTransformationMatrix3DNull::MutateTransformationMatrix3DNull()
    {
    }

    //! @brief Clone is the virtual copy constructor
    //! @return new instance of this class
    MutateTransformationMatrix3DNull *MutateTransformationMatrix3DNull::Clone() const
    {
      return new MutateTransformationMatrix3DNull( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateTransformationMatrix3DNull::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator which mutates a math::TransformationMatrix3D - in this case it does nothing
    //! @param MATRIX the math::TransformationMatrix3D which is to be mutated
    //! @return returns a math::MutateResult with the mutated math::TransformationMatrix3D - in this case unchanged
    math::MutateResult< math::TransformationMatrix3D> MutateTransformationMatrix3DNull::operator()
    (
      const math::TransformationMatrix3D &MATRIX
    ) const
    {
      BCL_MessageDbg( "MutateTransformationMatrix3DNull::operator()\n");

      // util::ShPtr< math::TransformationMatrix3D> sp_matrix( MATRIX.Clone()); //< this does not work

      // create ShPtr to math::TransformationMatrix3D initialized with "MATRIX"
      util::ShPtr< math::TransformationMatrix3D> sp_matrix( new math::TransformationMatrix3D( MATRIX));

      // return a math::MutateResult containing the mutated (unchanged) "sp_matrix"
      return math::MutateResult< math::TransformationMatrix3D>( sp_matrix, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateTransformationMatrix3DNull::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateTransformationMatrix3DNull::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
