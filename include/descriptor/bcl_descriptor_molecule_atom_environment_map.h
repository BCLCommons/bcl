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

#ifndef BCL_DESCRIPTOR_MOLECULE_ATOM_ENVIRONMENT_MAP_H_
#define BCL_DESCRIPTOR_MOLECULE_ATOM_ENVIRONMENT_MAP_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_hash_map.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeAtomEnvironmentMap
    //! Generates an atom environment hashmap of molecule fragments and checks to see if all
    //! of the atom environments can be found in a reference hashmap.
    //! Returns 1 if all fragment atom environments can be referenced, 0 otherwise
    //!
    //! @see @link example_descriptor_molecule_atom_environment_map.cpp @endlink
    //! @author brownbp1
    //! @date Apr 12, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeAtomEnvironmentMap :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

      //! atom environments
      static storage::HashMap< std::string, size_t> s_AtomEnvironments;

      //! atom environments reference input filename
      std::string m_AtomEnvironmentsFilename;

      //! check atom environment strings sans hydrogen atoms
      bool m_RemoveH;

      //! number of bond lengths to extend local atom environments
      size_t m_EnvironmentSize;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      MoleculeAtomEnvironmentMap();

      //! @brief constructor
      MoleculeAtomEnvironmentMap
      (
        const std::string &ENVIRONMENT_FILE,
        const bool &REMOVE_H,
        const size_t &ENVIRONMENT_SIZE
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeAtomEnvironmentMap
      MoleculeAtomEnvironmentMap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        linal::VectorReference< float> &STORAGE
      );

      //! @brief get a mutex
//      static sched::Mutex& GetMutex();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class MoleculeAtomEnvironmentMap

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_ATOM_ENVIRONMENT_MAP_H_
