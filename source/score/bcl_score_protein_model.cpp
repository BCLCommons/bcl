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
#include "score/bcl_score_protein_model.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  /////////////////
  // data access //
  /////////////////

    //! @brief conversion to a string from a Type
    //! @param TYPE the type to get a string for
    //! @return a string representing that type
    const std::string &ProteinModel::GetTypeName( const Type &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Score (Structure)",
        "Score (Sequence|Structure)",
        "Score (Misc.)",
        "undefined",
        GetStaticClassName< ProteinModel>()
      };
      return s_descriptors[ size_t( TYPE)];
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
    //! @param PROTEIN_MODEL_LHS first protein model score class
    //! @param PROTEIN_MODEL_RHS second protein model score class
    //! @return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
    bool ProteinModelLessThan::operator()
    (
      const ProteinModel &PROTEIN_MODEL_LHS, const ProteinModel &PROTEIN_MODEL_RHS
    ) const
    {
      // if the type of LHS is less than type of RHS
      if( PROTEIN_MODEL_LHS.GetType() < PROTEIN_MODEL_RHS.GetType())
      {
        return true;
      }

      // if the types are the same
      if( PROTEIN_MODEL_LHS.GetType() == PROTEIN_MODEL_RHS.GetType())
      {
        // check readable scheme
        if( PROTEIN_MODEL_LHS.GetReadableScheme() < PROTEIN_MODEL_RHS.GetReadableScheme())
        {
          return true;
        }
      }

      // otherwise return false
      return false;
    }

    //! @brief return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
    //! @param PROTEIN_MODEL_LHS first protein model score class
    //! @param PROTEIN_MODEL_RHS second protein model score class
    //! @return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
    bool ProteinModelLessThan::operator()
    (
      const util::PtrInterface< ProteinModel> &PROTEIN_MODEL_LHS,
      const util::PtrInterface< ProteinModel> &PROTEIN_MODEL_RHS
    ) const
    {
      return operator ()( *PROTEIN_MODEL_LHS, *PROTEIN_MODEL_RHS);
    }

    //! @brief return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
    //! @param PROTEIN_MODEL_LHS first protein model score class
    //! @param PROTEIN_MODEL_RHS second protein model score class
    //! @return true if PROTEIN_MODEL_LHS is less than PROTEIN_MODEL_RHS
    bool ProteinModelLessThan::operator()
    (
      const util::PtrInterface< const ProteinModel> &PROTEIN_MODEL_LHS,
      const util::PtrInterface< const ProteinModel> &PROTEIN_MODEL_RHS
    ) const
    {
      return operator ()( *PROTEIN_MODEL_LHS, *PROTEIN_MODEL_RHS);
    }

  } // namespace score
} // namespace bcl
