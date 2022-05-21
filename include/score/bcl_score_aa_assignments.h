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

#ifndef BCL_SCORE_AA_ASSIGNMENTS_H_
#define BCL_SCORE_AA_ASSIGNMENTS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "function/bcl_function_binary_interface.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAAssignments
    //! @brief enumerate all pairwise aa assignment scores
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_score_aa_assignments.cpp @endlink
    //! @author woetzen
    //! @date Apr 10, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAAssignments :
      public util::Enumerate< util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >, AAAssignments>
    {
      friend class util::Enumerate< util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >, AAAssignments>;

    public:

    //////////
    // data //
    //////////

      const AAAssignment e_IDENTITY;
      const AAAssignment e_PAM100;
      const AAAssignment e_PAM120;
      const AAAssignment e_PAM160;
      const AAAssignment e_PAM250;
      const AAAssignment e_BLOSUM90;
      const AAAssignment e_BLOSUM80;
      const AAAssignment e_BLOSUM62;
      const AAAssignment e_BLOSUM45;
      const AAAssignment e_PHAT85;
      const AAAssignment e_PHAT80;
      const AAAssignment e_PHAT75;
      const AAAssignment e_PHAT70;
      const AAAssignment e_BLAST;
      const AAAssignment e_PSIPRED;
      const AAAssignment e_JUFO;
      const AAAssignment e_SAM;
      const AAAssignment e_TMHMM;
      const AAAssignment e_TMMOD;
      const AAAssignment e_B2TMPRED;
      const AAAssignment e_PROFTMB;
      const AAAssignment e_CONPRED;
      const AAAssignment e_STERICAL_PARAMETER;
      const AAAssignment e_POLARIZABILITY;
      const AAAssignment e_VOLUME;
      const AAAssignment e_HYDROPHOBICITY;
      const AAAssignment e_ISOELECTRIC_POINT;
      const AAAssignment e_TFE_WHITE;
      const AAAssignment e_TFE_ENGELMAN;

    private:

      //! map of z-scores
      storage::Map< AAAssignment, util::ShPtr< math::ZScore> > m_ZScores;

      //! map of scores that have additional files associated with them
      storage::Map< AAAssignment, std::string> m_ScoreFile;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAAssignments();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief z-score associated with that enum
      //! @param ASSIGNMENT_ENUM score
      //! @return z-score to normalize the result of the energy function
      const util::ShPtr< math::ZScore> &GetZScore( const AAAssignment &ASSIGNMENT_ENUM) const;

      //! @brief score file extension, if score requires file to be read
      //! @param ASSIGNMENT_ENUM score
      //! @return the file extension required for the score - empty if non is required
      const std::string &GetFileExtension( const AAAssignment &ASSIGNMENT_ENUM) const;

    ////////////////
    // operations //
    ////////////////

    }; // class AAAssignments

    //! @brief construct on access function for all AAAssignments
    //! @return reference to only instances of AAAssignments
    BCL_API
    AAAssignments &GetAAAssignments();

  } // namespace score

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >, score::AAAssignments>;

  } // namespace util
} // namespace bcl

#endif // BCL_SCORE_AA_ASSIGNMENTS_H_
