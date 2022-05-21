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

#ifndef BCL_IO_RETRIEVE_INTERFACE_H_
#define BCL_IO_RETRIEVE_INTERFACE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RetrieveInterface
    //! @brief class defines an interface for storing and retrieving objects of type t_DataType from data sources.
    //!
    //! @tparam t_DataType types can be eg. t_DataType objects (protein, smallmolecule), model::InterfacesMetaData
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Jan 09, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Type, typename t_EnsembleType>
    class RetrieveInterface :
      public virtual util::SerializableInterface
    {

    public:

      //! @brief Clone the RetrieveInterface
      //! @return pointer to new RetrieveInterface
      virtual RetrieveInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief number of t_Types in storage
      //! @return number of t_Types in this storage
      virtual size_t GetSize() const = 0;

      //! @brief get all keys for this storage
      //! @return all keys for this storage
      virtual storage::Vector< std::string> GetAllKeys() const = 0;

      //! @brief get the size of the given key (units defined by t_Type)
      //! @param KEY the identifier for the specific object
      //! @return the number of t_Type-defined units possessed by the key
      virtual size_t GetKeySize( const std::string &KEY) const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief get stored t_Type from key
      //! @param KEY key identifier for specific object
      //! @return t_Type of interest
      virtual t_Type Retrieve( const std::string &KEY) const = 0;

      //! @brief get molecule ensemble for given keys
      //! @param KEYS vector of keys
      //! @return t_EnsembleType objects corresponding to KEYS
      virtual t_EnsembleType RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const = 0;

      //! @brief get ensemble of stored molecules objects for a range of keys, specified by plain numbers
      //! @param RANGE range of keys
      //! @return t_EnsembleType corresponding to keys in RANGE
      virtual t_EnsembleType RetrieveEnsemble( const math::Range< size_t> &RANGE) const = 0;

    }; // class RetrieveInterface

  } // namespace io
} // namespace bcl

#endif // BCL_IO_RETRIEVE_INTERFACE_H_
