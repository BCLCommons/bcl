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

#ifndef BCL_DESCRIPTOR_UMOL2D_H_
#define BCL_DESCRIPTOR_UMOL2D_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_atom_environment_bender.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UMol2D
    //! @brief calculates the # of certain types of atom environments in the molecule
    //! @details UMol2D definition is from:
    //!   Veber, et. al. 2002, Molecular Properties That Influence the Oral Bioavailability of Drug Candidates
    //!
    //!   4 types of atom hashsing types
    //!   Element
    //!   Atom
    //!
    //! @see @link example_descriptor_umol2d.cpp @endlink
    //! @author vuot
    //! @date Jan 12, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API UMol2D :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {
    private:
      //typedef chemistry::AtomEnvironmentBender::Atom_type t_AtomType;
      typedef chemistry::AtomEnvironmentBender::AtomTypeEnum t_AtomTypeEnum;
      //! Single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      t_AtomTypeEnum m_AtomType;
      size_t m_AeNumber;
      size_t m_Sphere;
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      UMol2D();

      //! @brief Clone function
      //! @return pointer to new MoleculeRotatableBonds
      UMol2D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief return the name of the property without any parameters
      //! @return name of the property as string
      static storage::Vector< chemistry::AtomEnvironmentBender> GetAEs( const t_AtomTypeEnum &ATOM_TYPE, const size_t &SPHERE);

    }; // class Molprint2D

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_UMOL2D_H_

