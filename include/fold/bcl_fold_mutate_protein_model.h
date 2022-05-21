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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModel
    //! @brief class applies the provided move to the entire protein model
    //!
    //! @see @link example_fold_mutate_protein_model.cpp @endlink
    //! @author woetzen, karakam
    //! @date May 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModel :
      public math::MutateInterface< assemble::ProteinModel>
    {
    private:

    //////////
    // data //
    //////////

      //! moves the protein
      util::Implementation< coord::MoveInterface> m_Move;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModel();

      //! @brief constructor from a Move and a scheme
      //! @param MOVE function that performs the move on the protein
      //! @param SCHEME scheme to be used
      MutateProteinModel
      (
        const coord::MoveInterface &MOVE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModel>()
      );

      //! @brief clone
      MutateProteinModel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief = operator
      //! @param MUTATE MutateProteinModel to be assigned to
      //! @return this MutateProteinModel after being assigned to MUTATE
      MutateProteinModel &operator =( const MutateProteinModel &MUTATE);

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateProteinModel

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_PROTEIN_MODEL_H_
