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

#ifndef BCL_FOLD_MUTATE_MULTIMER_H_
#define BCL_FOLD_MUTATE_MULTIMER_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateMultimer
    //! @brief mutates a multimer
    //! @details translates the associated protein toward or away from the axis of symmetry to change the multimer
    //!
    //! @see @link example_fold_mutate_multimer.cpp @endlink
    //! @author weinerbe
    //! @date Nov 30, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateMultimer :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! whether symmetry is dihedral
      bool m_IsDihedral;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from members
      //! @param IS_DIHEDRAL is dihedral symmetry
      //! @param SCHEME scheme for the mutate
      MutateMultimer( const bool IS_DIHEDRAL = false, const std::string &SCHEME = GetStaticClassName< MutateMultimer>());

      //! @brief Clone function
      //! @return pointer to new MutateMultimer
      MutateMultimer *Clone() const;

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

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief performs the mutate on the multiplier associated with the passed protein model
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    }; // class MutateMultimer

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_MULTIMER_H_
