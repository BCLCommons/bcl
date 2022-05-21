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

#ifndef BCL_BIOL_RAMACHANDRAN_H_
#define BCL_BIOL_RAMACHANDRAN_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "random/bcl_random_histogram_2d_distribution.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Ramachandran
    //! @brief class for storing the Ramachandran phi/psi distributions and various convenience functions
    //! @details This class stores the Ramachandran plots that include the phi/psi distributions for each amino acid
    //! and secondary structure type. It also allows getting a random phi/psi pair for a given amino acid and sstype.
    //!
    //! @see @link example_biol_ramachandran.cpp @endlink
    //! @author karakam
    //! @date Dec 30, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Ramachandran :
      public util::SerializableInterface
    {

    private:

      //! typedef for Histogram2D and random::Histogram2DDistribution pair
      typedef storage::Pair< math::Histogram2D, random::Histogram2DDistribution> HistogramDistributionPair;

    //////////
    // data //
    //////////

      //! filename for the sstype_histogram_filename;
      std::string m_HistogramFilename;

      //! map that stores the phi/psi distributions and random generators for each amino acid
      storage::Map< AAType, HistogramDistributionPair> m_AATypeMap;

      //! map of maps that stores the phi/psi distributions and random generators for each SSType and amino acid
      storage::Map< SSType, storage::Map< AAType, HistogramDistributionPair> > m_SSTypeMap;

      //! map of maps that stores the phi/psi distributions and random generators for each SSType and amino acid in membrane
      //! environment
      storage::Map< SSType, storage::Map< AAType, HistogramDistributionPair> > m_SSTypeMapMembrane;

    public:

      //! @brief get the static default instance of this class
      //! @return static default instance of this class
      static const Ramachandran &GetDefaultInstance();

      //! @brief get the default SSTypeHistogramFilename
      //! @return the default SSTypeHistogramFilename
      static const std::string &GetDefaultHistogramFilename();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Ramachandran();

      //! @brief constructor from a AAType and SSType histogram filename
      //! @param SS_TYPE_HISTOGRAM_FILENAME filename for the phi/psi histogram according to SSTypes and AATypes
      Ramachandran( const std::string &SS_TYPE_HISTOGRAM_FILENAME);

      //! @brief Clone function
      //! @return pointer to new Ramachandran
      Ramachandran *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the HistogramFilename
      //! @return the HistogramFilename
      const std::string &GetHistogramFilename() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief return a random phi psi for the given AAType
      //! @param AA_TYPE AAType of interest
      //! @return pair of random phi and psi values for the given AAType
      storage::VectorND< 2, double> GetRandomPhiPsi( const AAType &AA_TYPE) const;

      //! @brief return a random phi psi for the given AAType and SSType
      //! @param AA_TYPE AAType of interest
      //! @param SS_TYPE SSType of interest
      //! @return pair of random phi and psi values for the given AAType and SSType
      storage::VectorND< 2, double> GetRandomPhiPsi
      (
        const AAType &AA_TYPE,
        const SSType SS_TYPE
      ) const;

      //! @brief return a random phi psi for the given AAType and SSType
      //! @param AA_TYPE AAType of interest
      //! @param SS_TYPE SSType of interest
      //! @return pair of random phi and psi values for the given AAType and SSType
      storage::VectorND< 2, double> GetRandomPhiPsi
      (
        const AAType &AA_TYPE,
        const SSType SS_TYPE,
        const EnvironmentType &ENV_TYPE
      ) const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief function to read the ramachandran distributions and initialize the class members
      void Initialize();

    }; // class Ramachandran

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_RAMACHANDRAN_H_
