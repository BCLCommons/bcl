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

#ifndef BCL_FOLD_PHI_PSI_GENERATOR_RAMACHANDRAN_H_
#define BCL_FOLD_PHI_PSI_GENERATOR_RAMACHANDRAN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ramachandran.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PhiPsiGeneratorRamachandran
    //! @brief is for generating a random pair of phi and psi angles as biased by the
    //! ramachandran distribution.
    //!
    //! @see @link example_fold_phi_psi_generator_ramachandran.cpp @endlink
    //! @author alexanns
    //! @date Sep 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PhiPsiGeneratorRamachandran :
      public math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
    {

    private:

    //////////
    // data //
    //////////

      //! the ramachandran distributions for each amino acid type
      util::Implementation< biol::Ramachandran> m_Ramachandran;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief get default instance
      //! @return default instance
      static const math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > &
      GetDefaultInstance();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PhiPsiGeneratorRamachandran();

      //! @brief constructor from a ramachandran object
      //! @param RAMACHANDRAN instance of Ramachandran to be used
      PhiPsiGeneratorRamachandran( const biol::Ramachandran &RAMACHANDRAN);

      //! @brief Clone function
      //! @return pointer to new PhiPsiGeneratorRamachandran
      PhiPsiGeneratorRamachandran *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking MutationResidue and returning storage::VectorND< 2, double> with phi psi, respectively
      //! @param RESIDUE MutationResidue whose phi or psi will be changed
      //! @return phi or psi value, respectively, randomly selected based on a ramachandran distribution
      storage::VectorND< 2, double> operator()( const MutationResidue &RESIDUE) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PhiPsiGeneratorRamachandran

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PHI_PSI_GENERATOR_RAMACHANDRAN_H_ 
