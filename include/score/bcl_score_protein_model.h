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

#ifndef BCL_SCORE_PROTEIN_MODEL_H_
#define BCL_SCORE_PROTEIN_MODEL_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_ptr_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModel
    //! @brief Interface class that takes a protein model and returns a score
    //!
    //! @remarks example unnecessary
    //! @author weinerbe
    //! @date Jun 28, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModel :
      public math::FunctionInterfaceSerializable< assemble::ProteinModel, double>
    {

    public:

    ///////////
    // enums //
    ///////////

      //! enum for score type
      enum Type
      {
        e_Structure,   // score depends on structure (i.e. SSE pair)
        e_Sequence,    // score depends on sequence (i.e. AA environment)
        e_Misc,        // misc score
        e_Undefined,
        s_NumberTypes
      };

      //! @brief conversion to a string from a Type
      //! @param TYPE the type to get a string for
      //! @return a string representing that type
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief enum class wrapper for Type
      typedef util::WrapperEnum< Type, &GetTypeName, s_NumberTypes> TypeEnum;

      //! @brief virtual copy constructor
      virtual ProteinModel *Clone() const = 0;

      //! @brief detect quality-based scores. These are necessarily excluded from some algorithms
      virtual bool IsQualityScore() const
      {
        return false;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      virtual const std::string &GetReadableScheme() const
      {
        return GetScheme();
      }

      //! @brief get score type
      //! @return score type
      virtual Type GetType() const
      {
        return e_Undefined;
      }

    }; // class ProteinModel

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelLessThan
    //! @brief compares two ProteinModel objects
    //!
    //! @author weinerbe
    //! @date Jun 30, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelLessThan
    {
    public:

    ///////////////
    // operators //
    ///////////////

      //! @brief return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
      //! @param PROTEIN_MODEL_LHS first protein model score class
      //! @param PROTEIN_MODEL_RHS second protein model score class
      //! @return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
      bool operator()( const ProteinModel &PROTEIN_MODEL_LHS, const ProteinModel &PROTEIN_MODEL_RHS) const;

      //! @brief return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
      //! @param PROTEIN_MODEL_LHS first protein model score class
      //! @param PROTEIN_MODEL_RHS second protein model score class
      //! @return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
      bool operator()
      (
        const util::PtrInterface< ProteinModel> &PROTEIN_MODEL_LHS,
        const util::PtrInterface< ProteinModel> &PROTEIN_MODEL_RHS
      ) const;

      //! @brief return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
      //! @param PROTEIN_MODEL_LHS first protein model score class
      //! @param PROTEIN_MODEL_RHS second protein model score class
      //! @return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
      bool operator()
      (
        const util::PtrInterface< const ProteinModel> &PROTEIN_MODEL_LHS,
        const util::PtrInterface< const ProteinModel> &PROTEIN_MODEL_RHS
      ) const;

    }; // class ProteinModelLessThan

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_H_
