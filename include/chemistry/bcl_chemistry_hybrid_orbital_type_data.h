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

#ifndef BCL_CHEMISTRY_HYBRID_ORBITAL_TYPE_DATA_H_
#define BCL_CHEMISTRY_HYBRID_ORBITAL_TYPE_DATA_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atomic_orbital_types.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HybridOrbitalTypeData
    //! @brief to store a hybrid orbital definition
    //! @details this class abstracts the definition of a hybrid orbital to the orbitals that are hybridized
    //!
    //! @see @link example_chemistry_hybrid_orbital_data.cpp @endlink
    //! @author mueller, woetzen
    //! @date Aug 23, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HybridOrbitalTypeData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::Set< AtomicOrbitalTypesEnum> m_Orbitals; //!< orbitals that contribute to the hybrid

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
      HybridOrbitalTypeData();

      //! @brief construct from orbitals
      //! @param ORBITALS orbitals that are part of the hybrid
      HybridOrbitalTypeData
      (
        const storage::Set< AtomicOrbitalTypesEnum> &ORBITALS
      );

      //! @brief Clone function
      //! @return pointer to new HybridOrbitalTypeData
      HybridOrbitalTypeData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief access to the orbitals
      //! @return set of orbitals that make up this hybrid orbital
      const storage::Set< AtomicOrbitalTypesEnum> &GetOrbitals() const;

      //! @brief number fo possible sigma bonding partners
      //! @return the number of possible bonding partners bonded by sigma bonds
      size_t GetNumberOfPossibleSigmaBondingPartners() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

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
      //! @param INDENT number of indentations
      //! @return outputstream which was Readwritten to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class HybridOrbitalTypeData

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_HYBRID_ORBITAL_TYPE_DATA_H_ 
