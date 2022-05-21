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
#include "chemistry/bcl_chemistry_hybrid_orbital_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> HybridOrbitalTypeData::s_Instance
    (
      GetObjectInstances().AddInstance( new HybridOrbitalTypeData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HybridOrbitalTypeData::HybridOrbitalTypeData()
    {
    }

    //! @brief construct from orbitals
    //! @param ORBITALS orbitals that are part of the hybrid
    HybridOrbitalTypeData::HybridOrbitalTypeData
    (
      const storage::Set< AtomicOrbitalTypesEnum> &ORBITALS
    ) :
      m_Orbitals( ORBITALS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new HybridOrbitalTypeData
    HybridOrbitalTypeData *HybridOrbitalTypeData::Clone() const
    {
      return new HybridOrbitalTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &HybridOrbitalTypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the orbitals
    //! @return set of orbitals that make up this hybrid orbital
    const storage::Set< AtomicOrbitalTypesEnum> &HybridOrbitalTypeData::GetOrbitals() const
    {
      return m_Orbitals;
    }

    //! @brief number fo possible sigma bonding partners
    //! @return the number of possible bonding partners bonded by sigma bonds
    size_t HybridOrbitalTypeData::GetNumberOfPossibleSigmaBondingPartners() const
    {
      return m_Orbitals.GetSize();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HybridOrbitalTypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Orbitals, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was Readwritten to
    std::ostream &HybridOrbitalTypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Orbitals, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
