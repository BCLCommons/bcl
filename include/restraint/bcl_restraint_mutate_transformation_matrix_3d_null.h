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

#ifndef BCL_RESTRAINT_MUTATE_TRANSFORMATION_MATRIX_3D_NULL_H_
#define BCL_RESTRAINT_MUTATE_TRANSFORMATION_MATRIX_3D_NULL_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateTransformationMatrix3DNull
    //! @brief MutateTransformationMatrix3DNull is a class for mutating a transformation matrix that does not actually
    //! do a mutation
    //!
    //! @see @link example_restraint_mutate_transformation_matrix_3d_null.cpp @endlink
    //! @author alexanns
    //! @date 07/01/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateTransformationMatrix3DNull :
      public math::MutateInterface< math::TransformationMatrix3D>
    {
    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateTransformationMatrix3DNull();

      //! @brief Clone is the virtual copy constructor
      //! @return new instance of this class
      virtual MutateTransformationMatrix3DNull *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator which mutates a math::TransformationMatrix3D - in this case it does nothing
      //! @param MATRIX the math::TransformationMatrix3D which is to be mutated
      //! @return returns a math::MutateResult with the mutated math::TransformationMatrix3D - in this case unchanged
      math::MutateResult< math::TransformationMatrix3D> operator()
      (
        const math::TransformationMatrix3D &MATRIX
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateTransformationMatrix3DNull

  } // namespace restraint
} // namespace bcl

#endif //BCL_RESTRAINT_MUTATE_TRANSFORMATION_MATRIX_3D_NULL_H_
