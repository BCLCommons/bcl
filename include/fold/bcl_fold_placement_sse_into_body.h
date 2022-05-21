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

#ifndef BCL_FOLD_PLACEMENT_SSE_INTO_BODY_H_
#define BCL_FOLD_PLACEMENT_SSE_INTO_BODY_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "find/bcl_find.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_null.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementSSEIntoBody
    //! @brief class for placing an SSE into a ProteinModel as dictated by a restraint::Body
    //! @details This class uses the member restraint body and the orientation to determine where and how to place a given SSE
    //!
    //! @see @link example_fold_placement_sse_into_body.cpp @endlink
    //! @author alexanns
    //! @date June 30, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementSSEIntoBody :
      public PlacementInterface< assemble::SSE, assemble::ProteinModel>
    {

    private:

      //! ShPtr to a restraint:Body which is the restraint which provides the bodies an sse could be placed into
      util::ShPtr< restraint::Body> m_BodyRestraint;

      //! ShPtr to a PickCriteriaInterface which determines how a assemble::SSEGeometryInterface is picked out of "m_BodyRestraint"
      //! which is then used as the assemble::SSEGeometryInterface an SSE is placed into
      util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
        >
      > m_BodyPicker;

      //! ShPtr to a MutateInterface which determines how the SSE is placed into the body
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > m_Orientation;

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
      PlacementSSEIntoBody();

      //! @brief construct from the three member types
      //! @param RESTRAINT ShPtr to a restraint::Body which will be "m_BodyRestraint"
      //! @param BODY_PICKER ShPtr to a PickCriteriaInterface which will be "m_BodyPicker"
      //! @param ORIENTATION ShPtr to a FunctionInterface which will be "m_Orientation"
      PlacementSSEIntoBody
      (
        const util::ShPtr< restraint::Body> &RESTRAINT,
        const util::ShPtr
        <
          find::PickCriteriaInterface
          <
            util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
          >
        > &BODY_PICKER,
        util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > &ORIENTATION
      );

      //! @brief virtual copy constructor
      PlacementSSEIntoBody *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetBodyRestraint returns a const reference to "m_BodyRestraint"
      //! @return returns a const reference to "m_BodyRestraint"
      const util::ShPtr< restraint::Body> &GetBodyRestraint() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Place returns transformatin matrix for the placement of an SSE into a ProteinModel
      //! @param SSE SSE which is to be placed
      //! @param PROTEIN_MODEL into which the SSE will be placed
      //! @return the transformationmatrix3d to place SSE in PROTEIN_MODEL
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SSE,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! returns placement for the t_ObjectType without taking ProteinModel into consideration i.e. if it is empty
      //! @param SSE SSE which is to be placed
      //! @return the transformationmatrix3d to place SSE in PROTEIN_MODEL
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SSE
      ) const;

      //! @brief DetermineTransformationMatrix3DandBool figures out what the transformation matrix and bool should be
      //! @param RESTRAINT_BODY that was chosen to be the body into which the SSE will be placed
      //! @return the transformationmatrix3d to place SSE in PROTEIN_MODEL
      storage::Pair< math::TransformationMatrix3D, bool> DetermineTransformatinMatrix3DandBool
      (
        const util::ShPtr< assemble::SSEGeometryInterface> &RESTRAINT_BODY
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT indentation
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PlacementSSEIntoBody

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_PLACEMENT_SSE_INTO_BODY_H_
